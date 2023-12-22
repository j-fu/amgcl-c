add_compile_options(${OpenMP_CXX_FLAGS})
add_link_options(${OpenMP_CXX_FLAGS})

add_library(amgcl-c SHARED amgcl-c.cxx)
target_include_directories(amgcl-c
                           PUBLIC
                           "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/../include"
                           "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_include_directories(amgcl-c SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})
        

