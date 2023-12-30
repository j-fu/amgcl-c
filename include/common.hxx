#include <amgcl/backend/builtin.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/solver/runtime.hpp>
#include <amgcl/relaxation/as_preconditioner.hpp>
#include <amgcl/adapter/crs_tuple.hpp>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <istream>
#include <iostream>
#include <regex>
#include <complex.h>

extern "C"
{
#include "amgcl_c.h"
}


boost::property_tree::ptree boost_params(char *params);

template<typename Tv, typename Ti>
auto make_matrix_tuple(Ti n, Ti *ia, Ti *ja, Tv *a)
{
  Ti nnz=ia[n];
  auto ptr=amgcl::make_iterator_range(ia, ia + n+1);
  auto col=amgcl::make_iterator_range(ja, ja + nnz);
  auto val=amgcl::make_iterator_range(a, a + n);
  return std::make_tuple(n, ptr, col, val);
}

template <typename S, typename T> void destroy(S solver)
{  
  delete static_cast<T*>(solver.handle);
}

template<typename S, typename T, typename Tv, typename Ti> S create(Ti n,Ti *ia, Ti *ja, Tv *a, char *params)
{
  S solver;
  solver.handle=static_cast<void*>(new T(make_matrix_tuple(n,ia,ja,a), boost_params(params) ));
  return solver;
}

template <typename S, typename T, typename Tv> amgclcInfo solve(S _solver, Tv *sol, Tv *rhs)
{
  amgclcInfo info;
  auto solver = static_cast<T*>(_solver.handle);
  
  auto n = amgcl::backend::rows(solver->system_matrix());
  auto Sol=amgcl::make_iterator_range(sol, sol + n);
  auto Rhs=amgcl::make_iterator_range(rhs, rhs + n);
  
  std::tie(info.iters, info.residual) = (*solver)(Rhs,Sol);
  
  return info;
}

template <typename S, typename T,typename Tv> void apply(S _solver, Tv *sol, Tv *rhs)
{
  auto solver = static_cast<T*>(_solver.handle);
  
  auto n = amgcl::backend::rows(solver->system_matrix());
  auto Sol=amgcl::make_iterator_range(sol, sol + n);
  auto Rhs=amgcl::make_iterator_range(rhs, rhs + n);
  
  solver->apply(Rhs,Sol);
}


//
// AMG preconditioned Krylov solver
// See https://amgcl.readthedocs.io/en/latest/design.html?highlight=runtime#runtime-interface
//
template <typename Tv>
using AMGSolver=amgcl::make_solver<
  amgcl::amg<
    amgcl::backend::builtin<Tv>,
    amgcl::runtime::coarsening::wrapper,
    amgcl::runtime::relaxation::wrapper
    >,
  amgcl::runtime::solver::wrapper< amgcl::backend::builtin<Tv>>
  >;

const char *amgsolverparams();


//
// Relaxation preconditioned Krylov solver
//
template <typename Tv>
using RLXSolver=amgcl::make_solver<
  amgcl::relaxation::as_preconditioner<
    amgcl::backend::builtin<Tv>,
    amgcl::runtime::relaxation::wrapper
        >,
  amgcl::runtime::solver::wrapper<amgcl::backend::builtin<Tv>>
  > ;

const char *rlxsolverparams();

//
// AMG preconditioner
//
template <typename Tv>
using AMGPrecon = amgcl::amg<
    amgcl::backend::builtin<Tv>,
    amgcl::runtime::coarsening::wrapper,
    amgcl::runtime::relaxation::wrapper
  >;

const char *amgpreconparams();

//
// Relaxation preconditioner
//
template <typename Tv>
using RLXPrecon=amgcl::relaxation::as_preconditioner<
  amgcl::backend::builtin<Tv>,
  amgcl::runtime::relaxation::wrapper
  >;

const char *rlxpreconparams();
