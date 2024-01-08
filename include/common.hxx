#include <istream>
#include <iostream>
#include <regex>
#include <exception>
#include <cstring>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/solver/runtime.hpp>
#include <amgcl/relaxation/as_preconditioner.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/value_type/static_matrix.hpp>
#include <amgcl/adapter/block_matrix.hpp>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "amgcl_c.h"

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

template<typename S> S initialize_solver(S& solver,int blocksize)
{
  solver.error_state=0;
  solver.handle=0;
  solver.blocksize=blocksize;
  return solver;
}

template<typename I> I initialize_info(I& info)
{
  info.error_state=0;
  info.iters=0;
  info.residual=0.0;
  return info;
}


template<typename S> S set_error(S& solver, std::exception &e)
{
  solver.error_state=1;
  return solver;
}

template<typename S> S set_error(S& solver, int state)
{
  solver.error_state=state;
  return solver;
}


template<typename S, typename T, typename Tv, typename Ti> S create(Ti n,Ti *ia, Ti *ja, Tv *a, char *params)
{
  S solver;
  initialize_solver(solver,1);
  try
  {
    solver.handle=static_cast<void*>(new T(make_matrix_tuple(n,ia,ja,a), boost_params(params) ));
  }
  catch (std::exception &e)
  {
    set_error(solver,e);
    return solver;
  }
  return solver;
}

template<typename S, typename T, typename Tv, typename Ti, int N> S block_create(Ti n,Ti *ia, Ti *ja, Tv *a, char *params)
{
  S solver;
  initialize_solver(solver,N);
  try
  {
    auto A=make_matrix_tuple(n,ia,ja,a);
    auto Ab = amgcl::adapter::block_matrix<amgcl::static_matrix<Tv, N, N>>(A);
    solver.handle=static_cast<void*>(new T(Ab, boost_params(params) ));
  }
  catch (std::exception &e)
  {
    set_error(solver,e);
  }
  return solver;
}

template <typename S, typename T, typename Tv> amgclcInfo solve(S _solver, Tv *sol, Tv *rhs)
{
  amgclcInfo info;
  initialize_info(info);

  try
  {
    auto solver = static_cast<T*>(_solver.handle);
    auto n = amgcl::backend::rows(solver->system_matrix());
    auto Sol=amgcl::make_iterator_range(sol, sol + n);
    auto Rhs=amgcl::make_iterator_range(rhs, rhs + n);
    std::tie(info.iters, info.residual) = (*solver)(Rhs,Sol);
  }
  catch (std::exception &e)
  {
    set_error(info,e);
  }
  return info;
}


template <typename S, typename T, typename Tv, int N> amgclcInfo block_solve(S _solver, Tv *sol, Tv *rhs)
{
  amgclcInfo info;
  initialize_info(info);

  try
  {
    auto solver = static_cast<T*>(_solver.handle);
    auto bsol=reinterpret_cast<amgcl::static_matrix<double, N, 1>*>(sol);
    auto brhs=reinterpret_cast<amgcl::static_matrix<double, N, 1>*>(rhs);
    auto n = amgcl::backend::rows(solver->system_matrix());
    auto Sol=amgcl::make_iterator_range(bsol, bsol + n);
    auto Rhs=amgcl::make_iterator_range(brhs, brhs + n);
    
    std::tie(info.iters, info.residual) = (*solver)(Rhs,Sol);
  }
  catch (std::exception &e)
  {
    set_error(info,e);
  }
  
  return info;
}

template <typename S, typename T,typename Tv> amgclcInfo apply(S _solver, Tv *sol, Tv *rhs)
{
  amgclcInfo info;
  initialize_info(info);
  
  try
  {
    auto solver = static_cast<T*>(_solver.handle);
    auto n = amgcl::backend::rows(solver->system_matrix());
    auto Sol=amgcl::make_iterator_range(sol, sol + n);
    auto Rhs=amgcl::make_iterator_range(rhs, rhs + n);
    solver->apply(Rhs,Sol);
  }
  catch(std::exception &e)
  {
    set_error(info,e);
  }
  return info;
}


template <typename S, typename T,typename Tv, int N>  amgclcInfo block_apply(S _solver, Tv *sol, Tv *rhs)
{
  amgclcInfo info;
  initialize_info(info);
  
  try
  {
    auto solver = static_cast<T*>(_solver.handle);
    auto bsol=reinterpret_cast<amgcl::static_matrix<double, N, 1>*>(sol);
    auto brhs=reinterpret_cast<amgcl::static_matrix<double, N, 1>*>(rhs);
    auto n = amgcl::backend::rows(solver->system_matrix());
    auto Sol=amgcl::make_iterator_range(bsol, bsol + n);
    auto Rhs=amgcl::make_iterator_range(brhs, brhs + n);
    solver->apply(Rhs,Sol);
  }
  catch(std::exception &e)
  {
    set_error(info,e);
  }
  
  return info;
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

extern const char *amgsolverparams;


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

extern const char *rlxsolverparams;

//
// AMG preconditioner
//
template <typename Tv>
using AMGPrecon = amgcl::amg<
    amgcl::backend::builtin<Tv>,
    amgcl::runtime::coarsening::wrapper,
    amgcl::runtime::relaxation::wrapper
  >;

extern const char *amgpreconparams;

//
// Relaxation preconditioner
//
template <typename Tv>
using RLXPrecon=amgcl::relaxation::as_preconditioner<
  amgcl::backend::builtin<Tv>,
  amgcl::runtime::relaxation::wrapper
  >;

extern const char *rlxpreconparams;


//
// Block AMG preconditioned Krylov solver
// See https://amgcl.readthedocs.io/en/latest/design.html?highlight=runtime#runtime-interface
// and https://amgcl.readthedocs.io/en/latest/tutorial/Serena.html
//
template <typename Tv, int N>
using BlockAMGSolver=amgcl::make_solver<
  amgcl::amg<
    amgcl::backend::builtin< amgcl::static_matrix<Tv, N, N> >,
   amgcl::runtime::coarsening::wrapper,
    amgcl::runtime::relaxation::wrapper
    >,
  amgcl::runtime::solver::wrapper< amgcl::backend::builtin< amgcl::static_matrix<Tv, N, N> >>
  >;


//
// Block Relaxation preconditioned Krylov solver
//
template <typename Tv, int N>
using BlockRLXSolver=amgcl::make_solver<
  amgcl::relaxation::as_preconditioner<
    amgcl::backend::builtin< amgcl::static_matrix<Tv, N, N> >,
    amgcl::runtime::relaxation::wrapper
    >,
  amgcl::runtime::solver::wrapper<amgcl::backend::builtin< amgcl::static_matrix<Tv, N, N> >>
  > ;

//
// Block AMG preconditioner
//
template <typename Tv, int N>
using BlockAMGPrecon = amgcl::amg<
  amgcl::backend::builtin< amgcl::static_matrix<Tv, N, N> >,
    amgcl::runtime::coarsening::wrapper,
  amgcl::runtime::relaxation::wrapper
  >;


//
// Block Relaxation preconditioner
//
template <typename Tv, int N>
using BlockRLXPrecon=amgcl::relaxation::as_preconditioner<
  amgcl::backend::builtin< amgcl::static_matrix<Tv, N, N> >,
  amgcl::runtime::relaxation::wrapper
  >;

