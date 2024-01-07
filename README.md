AMGCL_C
========

Alternative C API for [AMGCL](https://github.com/ddemidov/amgcl).

The aim is the ability to access single level relaxation preconditioned Krylov methods as well as the possibility to perform single preconditioning steps (both for AMG and for relaxation) via a C interface. Solver parameters are passed as JSON strings.
It is also planned to build a Julia wrapper to AMGCL based on this code.

## API Description
Currently, AMGCL_C provides two interfaces:
- for `double` numbers and 4 byte int indexes, with the [hungarian notation](https://en.wikipedia.org/wiki/Hungarian_notation) prefix `amgclcDI` as shown in the example below
- for `double` numbers and 8 byte int indexes, with the [hungarian notation](https://en.wikipedia.org/wiki/Hungarian_notation) prefix `amgclcDL`
 just replace `int` with the coresponding 8 byte integer type (`long` or `long long` or `int64_t`).
  
### General parameters:
  - `...Create` parameters:
     -  `n,ia,ja,a`: Zero-based indexing CRS sparse matrix data
     -  `blocksize`: If blocksize >0, group unknowns into blocks of given size and cast the matrix internally to a sparse matrix of
        `blocksize x blocksize` static matrices. By default, block sizes 1...8 are instantiated. 
     -  `params`: JSON string containing parameter data. For easier escaping, `'` characters can be used as string delimiters instead of `"`. After the corresponding replacement, the string is parsed  by    [`boost::property_tree::json_parser::read_json`](https://www.boost.org/doc/libs/release/libs/property_tree/). This parser aborts upon syntax errors. When calling from C++, raw strings can be used to avoid the need of escaping end of line and `"` characters.
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
Please note that on Apple systems, OpenMP probably will not work (unless CMake
finds OpenMP).

```c
typedef struct{ void *handle; int blocksize;} amgclcDIAMGSolver;
amgclcDIAMGSolver amgclcDIAMGSolverCreate(int n, int *ia, int *ja, double *a, int blocksize, char *params);
amgclcInfo amgclcDIAMGSolverApply(amgclcDIAMGSolver solver, double *sol, double *rhs);
void amgclcDIAMGSolverDestroy(amgclcDIAMGSolver solver);
```

### Algebraic multigrid (AMG) preconditioned Krylov subspace iterative solver.

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
typedef struct{ void *handle; int blocksize;} amgclcDIRLXSolver;
amgclcDIRLXSolver amgclcDIRLXSolverCreate(int n, int *ia, int *ja, double *a, int blocksize, char *params);
amgclcInfo amgclcDIRLXSolverApply(amgclcDIRLXSolver solver, double *sol, double *rhs);
void amgclcDIRLXSolverDestroy(amgclcDIRLXSolver solver);
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
typedef struct{ void *handle; int blocksize;} amgclcDIAMGPrecon;
amgclcDIAMGPrecon amgclcDIAMGPreconCreate(int n, int *ia, int *ja, double *a, int blocksize, char *params);
void amgclcDIAMGPreconApply(amgclcDIAMGPrecon solver, double *sol, double *rhs);
void amgclcDIAMGPreconDestroy(amgclcDIAMGPrecon solver);
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
typedef struct{ void *handle; int blocksize;} amgclcDIRLXPrecon;
amgclcDIRLXPrecon amgclcDIRLXPreconCreate(int n, int *ia, int *ja, double *a, int blocksize, char *params);
void amgclcDIRLXPreconApply(amgclcDIRLXPrecon solver, double *sol, double *rhs);
void amgclcDIRLXPreconDestroy(amgclcDIRLXPrecon solver);
```

Default parameters:
```javascript
{
   "type": "ilu0"
}
```
