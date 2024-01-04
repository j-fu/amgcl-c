#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <amgcl_c.h>


#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>

#define TEST(expr) printf("TEST: %40s:",#expr), assert(expr), printf(" ok\n");

/* My own sparse matrix-vector multiplication for cross-checking amgcl claims */
void matmul(int n, int nnz, int *ia, int *ja, double *a, double *u, double *v)
{
  int i,j;
  for (i=0;i<n;i++)
  {
    v[i]=0;
    for (j=ia[i];j<ia[i+1];j++)
      v[i]+=a[j]*u[ja[j]];
  }
}

/* 3D poisson assembly, inspired from https://github.com/ddemidov/amgcl_benchmarks/blob/master/shared_mem/poisson/amgcl.cpp*/
int sample_problem(int n, int *_n3, int *_nnz, int **_ia, int **_ja, double **_a)
{
  int *ia, *ja;
  double *a;
  int i,j,k,idx;
  int cia,cja;
  
  int n3;
  double one,anisotropy, hx,hy,hz;
  one=1.0;
  anisotropy=1.0;
  n3=n*n*n;
  
  
  ia=(int*)calloc((n3+1),sizeof(int));
  ja=(int*)calloc((n3*7),sizeof(int));
  a=(double*)calloc((n3*7),sizeof(double));
  
  hx=1;
  hy=hx*anisotropy;
  hz=hy*anisotropy;

  cia=0;
  cja=0;

  ia[cia++]=0;
  for(k = 0, idx = 0; k < n; ++k)
  {
    for(j = 0; j < n; ++j)
    {
      for (i = 0; i < n; ++i, ++idx) {
        if (k > 0) {
          ja[cja]=(idx - n * n);
          a[cja]=(-1.0/(hz * hz) * one);
          cja++;
        }
        
        if (j > 0) {
          ja[cja]=(idx - n);
          a[cja]=(-1.0/(hy * hy) * one);
          cja++;
        }
        
        if (i > 0) {
          ja[cja]=(idx - 1);
          a[cja]=(-1.0/(hx * hx) * one);
          cja++;
        }
        
        ja[cja]=(idx);
        a[cja]=((2 / (hx * hx) + 2 / (hy * hy) + 2 / (hz * hz)) * one);
        cja++;
        if (i + 1 < n) {
          ja[cja]=(idx + 1);
          a[cja]=(-1.0/(hx * hx) * one);
          cja++;
        }

        if (j + 1 < n) {
          ja[cja]=(idx + n);
          a[cja]=(-1.0/(hy * hy) * one);
          cja++;
        }
        
        if (k + 1 < n) {
          ja[cja]=(idx + n * n);
          a[cja]=(-1.0/(hz * hz) * one);
          cja++;
        }
        
        ia[cia]=(cja);
        cia++;
      }
    }
  }
  *_n3=n3;
  *_nnz=cja;
  *_ia=ia;
  *_ja=ja;
  *_a=a;
}


int run(int n0, int blocksize)
{
  int n,nnz;
  int*ia,*ja;
  double *a,*rhs;
  double *u0, *u,*v;
  double myresidual,myresidual0;
  amgclcInfo info;
  amgclcDIAMGSolver amgsolver;
  amgclcDIRLXSolver rlxsolver;
  amgclcDIRLXPrecon rlxprecon;
  amgclcDIAMGPrecon amgprecon;
  
  int i;
  
  sample_problem(n0, &n, &nnz, &ia, &ja, &a);
  printf("n0=%d n=%d nnz=%d blocksize=%d\n",n0,n,nnz,blocksize);
  u=(double*)calloc(n,sizeof(double));
  u0=(double*)calloc(n,sizeof(double));
  v=(double*)calloc(n,sizeof(double));

  /*
    AMG solver
   */
  for(i=0;i<n;i++)
  {
    u[i]=1.0;
    v[i]=1.0;
  }
  amgsolver=amgclcDIAMGSolverCreate(n,ia,ja,a,blocksize,NULL);
  info=amgclcDIAMGSolverApply(amgsolver,u,v);
  amgclcDIAMGSolverDestroy(amgsolver);
  matmul(n,nnz,ia,ja,a,u,v);

  myresidual=0.0;
  for(i=0;i<n;i++)
  {
    myresidual+=(v[i]-1.0)*(v[i]-1.0)/n;
  }
  myresidual=sqrt(myresidual);
  printf("\namg solver: iters=%d residual=%e myresidual=%e\n",info.iters,info.residual,myresidual);
  if (n0==60)
  {
    TEST(info.residual<1.0e-10);
    TEST(myresidual<1.0e-10);
    TEST( fabs(myresidual-info.residual)<1.0e-14);
  }


  /*
    Relaxation solver
   */
  for(i=0;i<n;i++)
  {
    u[i]=1.0;
    v[i]=1.0;
  }
  rlxsolver=amgclcDIRLXSolverCreate(n,ia,ja,a,blocksize,NULL);
  info=amgclcDIRLXSolverApply(rlxsolver,u,v);
  amgclcDIRLXSolverDestroy(rlxsolver);
  matmul(n,nnz,ia,ja,a,u,v);

  myresidual=0.0;
  for(i=0;i<n;i++)
  {
    myresidual+=(v[i]-1.0)*(v[i]-1.0)/n;
  }
  myresidual=sqrt(myresidual);
  printf("\nrlx solver: iters=%d residual=%e myresidual=%e\n",info.iters,info.residual,myresidual);
  if (n0==60)
  {
    TEST(info.residual<1.0e-10);
    TEST(myresidual<1.0e-10);
    TEST( fabs(myresidual-info.residual)<1.0e-14);
  }

  /*
    AMG preconditioning step
   */
  for(i=0;i<n;i++)
  {
    u0[i]=1.0;
  }

  matmul(n,nnz,ia,ja,a,u0,v);
  myresidual0=0.0;
  for(i=0;i<n;i++)
  {
    myresidual0+=(v[i]-1.0)*(v[i]-1.0)/n;
    v[i]-=1.0;
  }
  myresidual0=sqrt(myresidual0);
  
  char *amgpreconparams="{'coarsening': { 'type': 'smoothed_aggregation'},   'relax': {'type': 'ilu0'} }";
  
  amgprecon=amgclcDIAMGPreconCreate(n,ia,ja,a,blocksize,amgpreconparams);
  amgclcDIAMGPreconApply(amgprecon,u,v);
  amgclcDIAMGPreconDestroy(amgprecon);
  for(i=0;i<n;i++)
    u[i]=u0[i]-u[i];
  
  matmul(n,nnz,ia,ja,a,u,v);
  myresidual=0.0;
  for(i=0;i<n;i++)
  {
    myresidual+=(v[i]-1.0)*(v[i]-1.0)/n;
  }
  myresidual=sqrt(myresidual);
  printf("\namg precon: myresidual/myresidual0=%e\n",myresidual/myresidual0);

  if (n0==60)
  {
    TEST(myresidual/myresidual0<0.8);
  }

  /*
    Relaxation precondtioning step for solving 
    Au=v with v[i]=1 and u0[i]=1
   */

  //  r= Au0-v
  for(i=0;i<n;i++)
  {
    u0[i]=1.0;
  }
  matmul(n,nnz,ia,ja,a,u0,v);
  myresidual0=0.0;
  for(i=0;i<n;i++)
  {
    myresidual0+=(v[i]-1.0)*(v[i]-1.0)/n;
    v[i]-=1.0;
  }
  myresidual0=sqrt(myresidual0);

  //m u=r
  rlxprecon=amgclcDIRLXPreconCreate(n,ia,ja,a,blocksize,NULL);
  amgclcDIRLXPreconApply(rlxprecon,u,v);
  
  amgclcDIRLXPreconDestroy(rlxprecon);
  for(i=0;i<n;i++)
    u[i]=u0[i]-u[i];

  matmul(n,nnz,ia,ja,a,u,v);
  myresidual=0.0;
  for(i=0;i<n;i++)
  {
    myresidual+=(v[i]-1.0)*(v[i]-1.0)/n;
  }
  myresidual=sqrt(myresidual);
  printf("\nrlx precon: myresidual/myresidual0=%e\n",myresidual/myresidual0);
  if (n0==60)
  {
    TEST(myresidual/myresidual0<1.0);
  }

  free(a);
  free(ia);
  free(ja);
  free(u);
  free(u0);
  free(v);
  return 0;
}

int main(int argc, char** argv)
{
  int n0;
  n0=60;
  if (argc>1)
  {
    n0=atoi(argv[1]);
  }
  run(n0,1);
  n0%2==0 && run(n0,2);
  n0%3==0 && run(n0,3);
  n0%4==0 && run(n0,4);
  n0%5==0 && run(n0,5);
  n0%6==0 && run(n0,6);
  n0%7==0 && run(n0,7);
  n0%8==0 && run(n0,8);

}

