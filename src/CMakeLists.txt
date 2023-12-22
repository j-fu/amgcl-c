add_compile_options(${OpenMP_CXX_FLAGS})
add_link_options(${OpenMP_CXX_FLAGS})

add_library(amgcl-c SHARED amgcl-c.cxx)
target_include_directories(amgcl-c SYSTEM PUBLIC ${AMGCL_DIR})
target_include_directories(amgcl-c SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})
        

