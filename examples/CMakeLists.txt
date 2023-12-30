if (OpenMP_FOUND)
  add_compile_options(${OpenMP_CXX_FLAGS})
  add_link_options(${OpenMP_CXX_FLAGS})
endif()

add_executable(example_capi example_capi.c)
add_dependencies(example_capi amgcl_c)

if (OpenMP_FOUND)
  target_link_libraries(example_capi amgcl_c m OpenMP::OpenMP_CXX)
else()
  target_link_libraries(example_capi amgcl_c m)
endif()

add_test(NAME example_capi COMMAND example_capi)

