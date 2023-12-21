#ifndef AMGCL_C_H
#define AMGCL_C_H

typedef struct {
  int iters;
  double error;
}  amgclcInfo;


typedef struct{
  void *handle;
} amgclcDAMGSolver;


amgclcDAMGSolver amgclcDAMGSolverCreate(int n, int *ia, int *ja, double *a, char *params);
amgclcInfo amgclcDAMGSolverApply(amgclcDAMGSolver solver, double *sol, double *rhs);
void amgclcDAMGSolverDestroy(amgclcDAMGSolver solver);

#endif
