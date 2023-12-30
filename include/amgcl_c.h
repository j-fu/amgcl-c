#ifndef AMGCL_C_H
#define AMGCL_C_H

#ifdef __cplusplus
extern "C" {
#endif

/**
  \file amgcl_c.h
  \brief Yet another C-API for amgcl_c.
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

/****************************************************************************/
/* DOUBLE/INT API */
  
/** 
    \brief algebraic multigrid preconditioned Krylov subspace iterative solver.

*/
typedef struct{ void *handle;} amgclcDIAMGSolver;
amgclcDIAMGSolver amgclcDIAMGSolverCreate(int n, int *ia, int *ja, double *a, char *params);
amgclcInfo amgclcDIAMGSolverApply(amgclcDIAMGSolver solver, double *sol, double *rhs);
void amgclcDIAMGSolverDestroy(amgclcDIAMGSolver solver);

/** 
    \brief Single level relaxation preconditioned Krylov subspace iterative solver.

 */
typedef struct{ void *handle;} amgclcDIRLXSolver;
amgclcDIRLXSolver amgclcDIRLXSolverCreate(int n, int *ia, int *ja, double *a, char *params);
amgclcInfo amgclcDIRLXSolverApply(amgclcDIRLXSolver solver, double *sol, double *rhs);
void amgclcDIRLXSolverDestroy(amgclcDIRLXSolver solver);

/** 
    \brief One algebraic multigrid preconditioning step.

 */
typedef struct{ void *handle;} amgclcDIAMGPrecon;
amgclcDIAMGPrecon amgclcDIAMGPreconCreate(int n, int *ia, int *ja, double *a, char *params);
void amgclcDIAMGPreconApply(amgclcDIAMGPrecon solver, double *sol, double *rhs);
void amgclcDIAMGPreconDestroy(amgclcDIAMGPrecon solver);

/** 
    \brief Ome single level relaxation  preconditioning step.

 */
typedef struct{ void *handle;} amgclcDIRLXPrecon;
amgclcDIRLXPrecon amgclcDIRLXPreconCreate(int n, int *ia, int *ja, double *a, char *params);
void amgclcDIRLXPreconApply(amgclcDIRLXPrecon solver, double *sol, double *rhs);
void amgclcDIRLXPreconDestroy(amgclcDIRLXPrecon solver);



/****************************************************************************/
/* DOUBLE/LONG API */

/** 
    \brief algebraic multigrid preconditioned Krylov subspace iterative solver.

*/
typedef struct{ void *handle;} amgclcDLAMGSolver;
amgclcDLAMGSolver amgclcDLAMGSolverCreate(long n, long *ia, long *ja, double *a, char *params);
amgclcInfo amgclcDLAMGSolverApply(amgclcDLAMGSolver solver, double *sol, double *rhs);
void amgclcDLAMGSolverDestroy(amgclcDLAMGSolver solver);

/** 
    \brief Single level relaxation preconditioned Krylov subspace iterative solver.

 */
typedef struct{ void *handle;} amgclcDLRLXSolver;
amgclcDLRLXSolver amgclcDLRLXSolverCreate(long n, long *ia, long *ja, double *a, char *params);
amgclcInfo amgclcDLRLXSolverApply(amgclcDLRLXSolver solver, double *sol, double *rhs);
void amgclcDLRLXSolverDestroy(amgclcDLRLXSolver solver);

/** 
    \brief One algebraic multigrid preconditioning step.

 */
typedef struct{ void *handle;} amgclcDLAMGPrecon;
amgclcDLAMGPrecon amgclcDLAMGPreconCreate(long n, long *ia, long *ja, double *a, char *params);
void amgclcDLAMGPreconApply(amgclcDLAMGPrecon solver, double *sol, double *rhs);
void amgclcDLAMGPreconDestroy(amgclcDLAMGPrecon solver);

/** 
    \brief Ome single level relaxation  preconditioning step.

 */
typedef struct{ void *handle;} amgclcDLRLXPrecon;
amgclcDLRLXPrecon amgclcDLRLXPreconCreate(long n, long *ia, long *ja, double *a, char *params);
void amgclcDLRLXPreconApply(amgclcDLRLXPrecon solver, double *sol, double *rhs);
void amgclcDLRLXPreconDestroy(amgclcDLRLXPrecon solver);


  
#ifdef __cplusplus
} // extern "C"
#endif

#endif
