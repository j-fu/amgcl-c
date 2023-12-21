#include <amgcl/backend/builtin.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/solver/runtime.hpp>
#include <amgcl/adapter/crs_tuple.hpp>

#include <boost/property_tree/json_parser.hpp>

#include <istream>
#include <iostream>
#include <regex>

extern "C"
{
#include "amgcl-c.h"
}


boost::property_tree::ptree boost_params(char *params)
{
  boost::property_tree::ptree prm;
  std::stringstream ssparams(std::regex_replace(std::string(params), std::regex("\'"), "\""));
  boost::property_tree::json_parser::read_json(ssparams,prm);

  // std::ostringstream os;
  // boost::property_tree::json_parser::write_json(os,prm);
  // std::cout << os.str() << std::endl;
  return prm;
}

auto make_matrix_tuple(int n,int *ia, int *ja, double *a)
{
  int nnz=ia[n+1];
  auto ptr=amgcl::make_iterator_range(ia, ia + n+1);
  auto col=amgcl::make_iterator_range(ja, ja + nnz);
  auto val=amgcl::make_iterator_range(a, a + n);
  return std::make_tuple(n, ptr, col, val);
}

template <typename S, typename T> void destroy(S solver)
{  
  delete static_cast<T*>(solver.handle);
}

template <typename S, typename T> S create(int n,int *ia, int *ja, double *a, char *params)
{
  S solver;
  solver.handle=static_cast<void*>(new T(make_matrix_tuple(n,ia,ja,a), boost_params(params) ));
  return solver;
}

template <typename S, typename T> amgclcInfo apply(S _solver, double *sol, double *rhs)
{
  amgclcInfo info;
  auto solver = static_cast<T*>(_solver.handle);
  
  auto n = amgcl::backend::rows(solver->system_matrix());
  auto Sol=amgcl::make_iterator_range(sol, sol + n);
  auto Rhs=amgcl::make_iterator_range(rhs, rhs + n);
  
  std::tie(info.iters, info.error) = (*solver)(Rhs,Sol);
  
  return info;
}


//
//  Built-in backend for doubles
// 
typedef amgcl::backend::builtin<double> DBuiltinBackend;


//
// AMG Solver
// See https://amgcl.readthedocs.io/en/latest/design.html?highlight=runtime#runtime-interface
//
typedef amgcl::make_solver<
  amgcl::amg<
    DBuiltinBackend,
    amgcl::runtime::coarsening::wrapper,
    amgcl::runtime::relaxation::wrapper
    >,
  amgcl::runtime::solver::wrapper<DBuiltinBackend>
  > DAMGSolver;

amgclcDAMGSolver amgclcDAMGSolverCreate(int n,int *ia, int *ja, double *a,char *params)
{
  return create<amgclcDAMGSolver,DAMGSolver>(n,ia,ja,a,params);
}

amgclcInfo amgclcDAMGSolverApply(amgclcDAMGSolver solver, double *sol, double *rhs)
{
  return apply<amgclcDAMGSolver,DAMGSolver>(solver,sol,rhs);
}

void amgclcDAMGSolverDestroy(amgclcDAMGSolver solver)
{
  destroy<amgclcDAMGSolver,DAMGSolver>(solver);
}

