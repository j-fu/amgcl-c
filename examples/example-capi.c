#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <amgcl-c.h>




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

int main(int argc, char** argv)
{
  int n0,n,nnz;
  int*ia,*ja;
  double *a,*rhs;
  double *u,*v;
  double myresidual;
  amgclcInfo info;
  amgclcDAMGSolver amgsolver;
  amgclcDRLXSolver rlxsolver;
  amgclcDRLXPrecon rlxprecon;
  amgclcDAMGPrecon amgprecon;
  
  int i;
  printf("main:\n");

  n0=10;
  if (argc>1)
  {
    n0=atoi(argv[1]);
  }
  
  sample_problem(n0, &n, &nnz, &ia, &ja, &a);
  printf("n0=%d n=%d nnz=%d\n",n0,n,nnz);
  u=(double*)calloc(n,sizeof(double));
  v=(double*)calloc(n,sizeof(double));

  for(i=0;i<n;i++)
  {
    u[i]=1.0;
    v[i]=1.0;
  }

  char *amgparams="{\
    'solver': {\
      'type': 'bicgstab',\
      'tol': 1.0e-10,\
      'maxiter': 10},\
    'precond': {\
      'coarsening': {\
        'type': 'smoothed_aggregation', 'relax': 1.0},\
        'relax': {\
        'type': 'spai0'}\
    }\
}";
  

  
  amgsolver=amgclcDAMGSolverCreate(n,ia,ja,a,amgparams);
  info=amgclcDAMGSolverApply(amgsolver,u,v);
  amgclcDAMGSolverDestroy(amgsolver);
  matmul(n,nnz,ia,ja,a,u,v);

  myresidual=0.0;
  for(i=0;i<n;i++)
  {
    myresidual+=(v[i]-1.0)*(v[i]-1.0)/n;
  }
  myresidual=sqrt(myresidual);
  printf("amg: iters=%d residual=%e myresidual=%e\n",info.iters,info.residual,myresidual);




    char *rlxparams="{\
    'solver': {'type': 'bicgstab','tol': 1.0e-10, 'maxiter': 100 },\
    'precond': {'type': 'ilu0' }\
    }";



  for(i=0;i<n;i++)
  {
    u[i]=1.0;
    v[i]=1.0;
  }

  rlxsolver=amgclcDRLXSolverCreate(n,ia,ja,a,rlxparams);
  info=amgclcDRLXSolverApply(rlxsolver,u,v);
  amgclcDRLXSolverDestroy(rlxsolver);
  matmul(n,nnz,ia,ja,a,u,v);

  myresidual=0.0;
  for(i=0;i<n;i++)
  {
    myresidual+=(v[i]-1.0)*(v[i]-1.0)/n;
  }
  myresidual=sqrt(myresidual);
  printf("rlx: iters=%d residual=%e myresidual=%e\n",info.iters,info.residual,myresidual);



  for(i=0;i<n;i++)
  {
    u[i]=1.0;
    v[i]=1.0;
  }

  rlxprecon=amgclcDRLXPreconCreate(n,ia,ja,a,rlxparams);
  amgclcDRLXPreconApply(rlxprecon,u,v);
  amgclcDRLXPreconDestroy(rlxprecon);
  matmul(n,nnz,ia,ja,a,u,v);

  myresidual=0.0;
  for(i=0;i<n;i++)
  {
    myresidual+=(v[i]-1.0)*(v[i]-1.0)/n;
  }
  myresidual=sqrt(myresidual);
  printf("rlxprecon: myresidual=%e\n",myresidual);


  for(i=0;i<n;i++)
  {
    u[i]=1.0;
    v[i]=1.0;
  }

  amgprecon=amgclcDAMGPreconCreate(n,ia,ja,a,amgparams);
  amgclcDAMGPreconApply(amgprecon,u,v);
  amgclcDAMGPreconDestroy(amgprecon);
  matmul(n,nnz,ia,ja,a,u,v);

  myresidual=0.0;
  for(i=0;i<n;i++)
  {
    myresidual+=(v[i]-1.0)*(v[i]-1.0)/n;
  }
  myresidual=sqrt(myresidual);
  printf("amgprecon: myresidual=%e\n",myresidual);



  
  free(a);
  free(ia);
  free(ja);
  free(u);
  free(v);
  printf("ok.\n");
  
}
