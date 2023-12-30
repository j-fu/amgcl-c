if (OpenMP_FOUND)
  add_compile_options(${OpenMP_CXX_FLAGS})
  add_link_options(${OpenMP_CXX_FLAGS})
endif()

add_library(amgcl_c SHARED amgcl_c.cxx)
target_include_directories(amgcl_c
                           PUBLIC
                           "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/../include"
                           "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_include_directories(amgcl_c SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})
        

