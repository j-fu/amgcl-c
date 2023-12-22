#ifndef AMGCL_C_H
#define AMGCL_C_H

#ifdef __cplusplus
extern "C" {
#endif

/**
  \file amgcl-c.h
  \brief Yet another C-API for amgcl-c.
  \author Jürgen Fuhrmann <juergen.fuhrmann@wias-berlin.de>
  
  - Create parameters:
     -  `n,ia,ja,a`: Zero-based indexing CRS sparse matrix data
     -  `params`: JSON string containing parameter data
     
  - Apply parameters:
     - sol, rhs: zero-offset vectors. Length is determined from the created 
       solver/preconditioner
*/

/**
  \brief Info structure returned by SolverApply functions.
*/
typedef struct {
  int iters;
  double residual;
}  amgclcInfo;


/** 
    \brief algebraic multigrid preconditioned Krylov subspace iterative solver.

*/
typedef struct{ void *handle;} amgclcDAMGSolver;
amgclcDAMGSolver amgclcDAMGSolverCreate(int n, int *ia, int *ja, double *a, char *params);
amgclcInfo amgclcDAMGSolverApply(amgclcDAMGSolver solver, double *sol, double *rhs);
void amgclcDAMGSolverDestroy(amgclcDAMGSolver solver);

/** 
    \brief Single level relaxation preconditioned Krylov subspace iterative solver.

 */
typedef struct{ void *handle;} amgclcDRLXSolver;
amgclcDRLXSolver amgclcDRLXSolverCreate(int n, int *ia, int *ja, double *a, char *params);
amgclcInfo amgclcDRLXSolverApply(amgclcDRLXSolver solver, double *sol, double *rhs);
void amgclcDRLXSolverDestroy(amgclcDRLXSolver solver);

/** 
    \brief One algebraic multigrid preconditioning step.

 */
typedef struct{ void *handle;} amgclcDAMGPrecon;
amgclcDAMGPrecon amgclcDAMGPreconCreate(int n, int *ia, int *ja, double *a, char *params);
void amgclcDAMGPreconApply(amgclcDAMGPrecon solver, double *sol, double *rhs);
void amgclcDAMGPreconDestroy(amgclcDAMGPrecon solver);

/** 
    \brief Ome single level relaxation  preconditioning step.

 */
typedef struct{ void *handle;} amgclcDRLXPrecon;
amgclcDRLXPrecon amgclcDRLXPreconCreate(int n, int *ia, int *ja, double *a, char *params);
void amgclcDRLXPreconApply(amgclcDRLXPrecon solver, double *sol, double *rhs);
void amgclcDRLXPreconDestroy(amgclcDRLXPrecon solver);

#ifdef __cplusplus
} // extern "C"
#endif

#endif
