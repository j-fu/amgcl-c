#include <common.hxx>

xxxSolver xxxCreate(int n)
{
  xxxSolver S;
  S.blocksize=n;
  S.handle=0;
  return S;
}

boost::property_tree::ptree boost_params(char *params)
{
  boost::property_tree::ptree prm;
  std::stringstream ssparams(std::regex_replace(std::string(params), std::regex("\'"), "\""));
  try
  {
    boost::property_tree::json_parser::read_json(ssparams,prm);
  }
  catch(const std::exception& e)
  {
    std::cerr << e.what() << std::endl;
    std::cerr << "Error parsing json string:" <<std::endl << ssparams.str() << std::endl;
    throw e;
  }
  // std::ostringstream os;
  // boost::property_tree::json_parser::write_json(os,prm);
  // std::cout << os.str() << std::endl;
  return prm;
}

const char *amgsolverparams= R"(
{"solver": { "type": "bicgstab",  "tol": 1.0e-10, "maxiter": 50},
    "precond": {
      "coarsening": { "type": "smoothed_aggregation", "relax": 1.0},
      "relax": {"type": "spai0"}
     }
}
)";


const char *rlxsolverparams= R"(
{
 "solver": {"type": "bicgstab","tol": 1.0e-10, "maxiter": 100 },
 "precond": {"type": "ilu0" }
}
)";

const char *amgpreconparams=R"(
{
   "coarsening": { "type": "smoothed_aggregation", "relax": 1.0},
   "relax": {"type": "spai0"}
}
)";

const char *rlxpreconparams=R"(
{
   "type": "ilu0"
}
)";


