#ifndef AMGCL_C_H
#define AMGCL_C_H

typedef struct {
  int iters;
  double error;
}  amgclcInfo;


typedef struct{ void *handle;} amgclcDAMGSolver;
amgclcDAMGSolver amgclcDAMGSolverCreate(int n, int *ia, int *ja, double *a, char *params);
amgclcInfo amgclcDAMGSolverApply(amgclcDAMGSolver solver, double *sol, double *rhs);
void amgclcDAMGSolverDestroy(amgclcDAMGSolver solver);


typedef struct{ void *handle;} amgclcDRLXSolver;
amgclcDRLXSolver amgclcDRLXSolverCreate(int n, int *ia, int *ja, double *a, char *params);
amgclcInfo amgclcDRLXSolverApply(amgclcDRLXSolver solver, double *sol, double *rhs);
void amgclcDRLXSolverDestroy(amgclcDRLXSolver solver);


typedef struct{ void *handle;} amgclcDAMGPrecon;
amgclcDAMGPrecon amgclcDAMGPreconCreate(int n, int *ia, int *ja, double *a, char *params);
void amgclcDAMGPreconApply(amgclcDAMGPrecon solver, double *sol, double *rhs);
void amgclcDAMGPreconDestroy(amgclcDAMGPrecon solver);



typedef struct{ void *handle;} amgclcDRLXPrecon;
amgclcDRLXPrecon amgclcDRLXPreconCreate(int n, int *ia, int *ja, double *a, char *params);
void amgclcDRLXPreconApply(amgclcDRLXPrecon solver, double *sol, double *rhs);
void amgclcDRLXPreconDestroy(amgclcDRLXPrecon solver);


#endif
