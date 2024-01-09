amgcl-c
========

Alternative C API for a subset of [AMGCL](https://github.com/ddemidov/amgcl).

The aim is the ability to access single level relaxation preconditioned Krylov methods as well as the possibility to perform single preconditioning steps (both for AMG and for relaxation) via a C interface. Solver parameters are passed as JSON strings.
It is also planned to build a Julia wrapper to AMGCL based on this code.

## API Description

### General parameters:
  - `...Create` parameters:
     -  `n,ia,ja,a`: Zero-based indexing CRS sparse matrix data
     -  `params`: JSON string containing parameter data. For easier escaping, `'` characters
     can be used as string delimiters instead of `"`. After the corresponding replacement, the string is parsed  by    [`boost::property_tree::json_parser::read_json`](https://www.boost.org/doc/libs/release/libs/property_tree/). This parser aborts upon syntax errors. When calling from C++, raw strings can be used to avoid the need of escaping end of line and `"` characters.
  - `...Apply` parameters:
     - `sol, rhs`: zero-offset vectors. Length is determined from the created solver/preconditioner
  -  Data structure returned by iterative methods:
```c
typedef struct {
  int iters;
  double residual;
}  amgclcInfo;
```

The only backend supported in the moment is AMGCL's default OpenMP parallel backend.

### Algebraic multigrid (AMG) preconditioned Krylov subspace iterative solver.

```c
typedef struct{ void *handle;} amgclcDAMGSolver;
amgclcDAMGSolver amgclcDAMGSolverCreate(int n, int *ia, int *ja, double *a, char *params);
amgclcInfo amgclcDAMGSolverApply(amgclcDAMGSolver solver, double *sol, double *rhs);
void amgclcDAMGSolverDestroy(amgclcDAMGSolver solver);
```

Default parameters:
```javascript
{
  "solver":  { "type": "bicgstab",  "tol": 1.0e-10, "maxiter": 10},
  "precond": {
       "coarsening": { "type": "smoothed_aggregation", "relax": 1.0},
        "relax":     {"type": "spai0"}
       }
}
```


### Single level relaxation preconditioned Krylov subspace iterative solver.

```c
typedef struct{ void *handle;} amgclcDRLXSolver;
amgclcDRLXSolver amgclcDRLXSolverCreate(int n, int *ia, int *ja, double *a, char *params);
amgclcInfo amgclcDRLXSolverApply(amgclcDRLXSolver solver, double *sol, double *rhs);
void amgclcDRLXSolverDestroy(amgclcDRLXSolver solver);
```

Default parameters:
```javascript
{
  "solver": {"type": "bicgstab","tol": 1.0e-10, "maxiter": 100 },
  "precond": {"type": "ilu0" }
}
```


### One AMG preconditioning step

```c
typedef struct{ void *handle;} amgclcDAMGPrecon;
amgclcDAMGPrecon amgclcDAMGPreconCreate(int n, int *ia, int *ja, double *a, char *params);
void amgclcDAMGPreconApply(amgclcDAMGPrecon solver, double *sol, double *rhs);
void amgclcDAMGPreconDestroy(amgclcDAMGPrecon solver);
```

Default parameters:

```javascript
{
   "coarsening": { "type": "smoothed_aggregation", "relax": 1.0},
   "relax": {"type": "spai0"}
}
```

### One single level relaxation  preconditioning step.

```c
typedef struct{ void *handle;} amgclcDRLXPrecon;
amgclcDRLXPrecon amgclcDRLXPreconCreate(int n, int *ia, int *ja, double *a, char *params);
void amgclcDRLXPreconApply(amgclcDRLXPrecon solver, double *sol, double *rhs);
void amgclcDRLXPreconDestroy(amgclcDRLXPrecon solver);
```

Default parameters:
```javascript
{
   "type": "ilu0"
}
```
