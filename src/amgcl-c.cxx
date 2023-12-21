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
typedef amgcl::backend::builtin<double> Backend;

typedef amgcl::make_solver<
    amgcl::amg<
        Backend,
        amgcl::runtime::coarsening::wrapper,
        amgcl::runtime::relaxation::wrapper
        >,
    amgcl::runtime::solver::wrapper<Backend>
    > Solver;


extern "C" amgclcDAMGSolver amgclcDAMGSolverCreate(int n,int *ia, int *ja, double *a,char *params)
{
  amgclcDAMGSolver solver;
  int nnz;
  
  nnz=ia[n+1];
  auto ptr=amgcl::make_iterator_range(ia, ia + n+1);
  auto col=amgcl::make_iterator_range(ja, ja + nnz);
  auto val=amgcl::make_iterator_range(a, a + n);
  
  boost::property_tree::ptree prm;
  
 
  std::stringstream ssparams(std::regex_replace(std::string(params), std::regex("\'"), "\""));
  
  boost::property_tree::json_parser::read_json(ssparams,prm);
  
  solver.handle=static_cast<void*>(new Solver( std::make_tuple(n, ptr, col, val), prm ));
  return solver;
}


extern "C" amgclcInfo amgclcDAMGSolverApply(amgclcDAMGSolver _solver, double *u, double *v)
{
  amgclcInfo info;
  Solver *solver = static_cast<Solver*>(_solver.handle);
  
  size_t n = amgcl::backend::rows(solver->system_matrix());
  auto U=amgcl::make_iterator_range(u, u + n);
  auto V=amgcl::make_iterator_range(v, v + n);
  
  std::tie(info.iters, info.error) = (*solver)(V, U);
  
  return info;
}


extern "C" void amgclcDAMGSolverDestroy(amgclcDAMGSolver solver)
{
  delete static_cast<Solver*>(solver.handle);
}
