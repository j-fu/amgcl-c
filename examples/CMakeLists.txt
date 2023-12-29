add_compile_options(${OpenMP_C_FLAGS})
add_link_options(${OpenMP_C_FLAGS})

add_executable(example_capi example_capi.c)
add_dependencies(example_capi amgcl_c)
target_link_libraries(example_capi amgcl_c m OpenMP::OpenMP_CXX)

add_test(NAME example_capi COMMAND example_capi)

