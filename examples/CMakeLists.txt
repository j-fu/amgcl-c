add_executable(example-capi example-capi.c)
add_dependencies(example-capi amgcl-c)
target_link_libraries(example-capi amgcl-c m)
