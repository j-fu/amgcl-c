#ifndef AMGCL_C_H
#define AMGCL_C_H

typedef struct {
  int iters;
  double error;
}  amgclc_info;


typedef void* amgclc_handle;

amgclc_handle create_solver(int n, int *ia, int *ja, double *a, char *params);
amgclc_info   apply_solver(amgclc_handle solver, double *u, double *v);
void destroy_solver(amgclc_handle solver);

#endif
