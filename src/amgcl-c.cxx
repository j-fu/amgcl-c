#include <amgcl/backend/builtin.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/solver/runtime.hpp>
#include <amgcl/adapter/crs_tuple.hpp>

typedef amgcl::backend::builtin<double> Backend;

typedef amgcl::make_solver<
    amgcl::amg<
        Backend,
        amgcl::runtime::coarsening::wrapper,
        amgcl::runtime::relaxation::wrapper
        >,
    amgcl::runtime::solver::wrapper<Backend>
    > Solver;


extern "C" void solve(int n, int nnz, int *ia, int *ja, double *a, double *u, double *v, int*iters, double *error)
{

  auto ptr=amgcl::make_iterator_range(ia, ia + n+1);
  auto col=amgcl::make_iterator_range(ja, ja + nnz);
  auto val=amgcl::make_iterator_range(a, a + n);
  auto U=amgcl::make_iterator_range(u, u + n);
  auto V=amgcl::make_iterator_range(v, v + n);

  
  boost::property_tree::ptree prm;
  
  prm.put("solver.type", "bicgstab");
  prm.put("solver.tol", 1e-3);
  prm.put("solver.maxiter", 10);
  prm.put("precond.coarsening.type", "smoothed_aggregation");
  prm.put("precond.relax.type", "spai0");

  Solver do_solve( std::make_tuple(n, ptr, col, val), prm );
  std::tie(*iters, *error) = do_solve(V, U);

}

