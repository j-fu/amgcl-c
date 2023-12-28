add_compile_options(${OpenMP_C_FLAGS})
add_link_options(${OpenMP_C_FLAGS})

add_executable(example-capi example-capi.c)
add_dependencies(example-capi amgcl-c)
target_link_libraries(example-capi amgcl-c m OpenMP::OpenMP_CXX)

add_test(NAME example-capi COMMAND example-capi)

