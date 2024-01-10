AMGCL_C
========

Alternative C API for a subset of [AMGCL](https://github.com/ddemidov/amgcl).

The aim is the ability to access AMG and single level relaxation preconditioned Krylov methods with the additional possibility to perform single preconditioning steps (both for AMG and for relaxation) via a C interface. Solver parameters are passed as JSON strings.

If you find this package useful, please  [give credits](https://github.com/ddemidov/amgcl?tab=readme-ov-file#referencing) to the author of AMGCL.

Compiled libraries for this package are available via the [AMGCL_C_jll.jl](https://github.com/JuliaBinaryWrappers/AMGCL_C_jll.jl) Julia  binary package which serves as the backend for the [AMGCLWrap.jl](https://github.com/j-fu/AMGCLWrap.jl) Julia package.


## API Description

For general functionality and possible parameter values, please consult the [AMGCL documentation](https://amgcl.readthedocs.io/en/latest/components.html).

Currently, AMGCL_C provides two interfaces:
- for `double` numbers and 4 byte int indexes, with the [hungarian notation](https://en.wikipedia.org/wiki/Hungarian_notation) prefix `amgclcDI` as shown in the example below
- for `double` numbers and 8 byte int indexes, with the [hungarian notation](https://en.wikipedia.org/wiki/Hungarian_notation) prefix `amgclcDL`.
  In the examples below just replace `DI` with `DL` and `int*` with a corresponding 8 byte pointer to integer type (`long*`,  `long long*` or `int64_t*`).
  
### General parameters:
  - `...Create` parameters:
     -  `n,ia,ja,a`: Zero-based indexing CRS sparse matrix data
     -  `blocksize`: If blocksize >0, group unknowns into blocks of given size and cast the matrix internally to a sparse matrix of `blocksize x blocksize` static matrices. By default, block sizes 1...8 are instantiated. Call `int amgclcBlocksizeInstantiated(blocksize)` to check if a given blocksize is instantiated. For uninstantiated block sizes, the `...Create` method returns with `error_state!=0`.
     -  `params`: JSON string containing parameter data. For easier escaping, `'` characters can be used as string delimiters instead of `"`. After the corresponding replacement, the string is parsed  by    [`boost::property_tree::json_parser::read_json`](https://www.boost.org/doc/libs/release/libs/property_tree/) into
a `boost::property_tree::p_tree` object which is fed into the [runtime interface](https://amgcl.readthedocs.io/en/latest/design.html?highlight=runtime#runtime-interface) of AMGCL.  Upon syntax errors, the runtime error thrown by the parser is detected, and the `...Create` method returns with `error_state!=0`. When calling from C++, raw strings can be used to avoid the need of escaping end of line and `"` characters.
  - `...Apply` parameters:
     - `sol, rhs`: zero-offset vectors. Length is determined from the created solver/preconditioner
  -  Data structure returned by iterative methods:
```c
typedef struct {
  int iters;
  double residual;
  int error_state;
}  amgclcInfo;
```

The `error_state` field in `amgclInfo` and the solver structs can be used to check for errors.

The only backend supported in the moment is AMGCL's default OpenMP parallel backend.

```c
typedef struct{ void *handle; int blocksize;int error_state;} amgclcDIAMGSolver;
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
typedef struct{ void *handle; int blocksize;int error_state;} amgclcDIRLXSolver;
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
typedef struct{ void *handle; int blocksize;int error_state;} amgclcDIAMGPrecon;
amgclcDIAMGPrecon amgclcDIAMGPreconCreate(int n, int *ia, int *ja, double *a, int blocksize, char *params);
amgclcInfo amgclcDIAMGPreconApply(amgclcDIAMGPrecon solver, double *sol, double *rhs);
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
typedef struct{ void *handle; int blocksize;int error_state;} amgclcDIRLXPrecon;
amgclcDIRLXPrecon amgclcDIRLXPreconCreate(int n, int *ia, int *ja, double *a, int blocksize, char *params);
amgclcInfo amgclcDIRLXPreconApply(amgclcDIRLXPrecon solver, double *sol, double *rhs);
void amgclcDIRLXPreconDestroy(amgclcDIRLXPrecon solver);
```

Default parameters:
```javascript
{
   "type": "ilu0"
}
```
