if (OpenMP_FOUND)
  add_compile_options(${OpenMP_CXX_FLAGS})
  add_link_options(${OpenMP_CXX_FLAGS})
endif()

add_custom_command(
    OUTPUT  amgcl_di.cxx
    COMMAND ${CMAKE_COMMAND} -D inputfile=${CMAKE_CURRENT_SOURCE_DIR}/amgcl_tv_ti.cxx.in -D BLOCKSIZES="${BLOCKSIZES}" -D outputfile=${CMAKE_CURRENT_BINARY_DIR}/amgcl_di.cxx  -D TVTI="DI"  -D TV="double"  -D TI="int"  -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/specialize.cmake
    DEPENDS amgcl_tv_ti.cxx.in ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/specialize.cmake amgcl_c_di.h
)

add_custom_command(
    OUTPUT  amgcl_dl.cxx
    COMMAND ${CMAKE_COMMAND} -D inputfile=${CMAKE_CURRENT_SOURCE_DIR}/amgcl_tv_ti.cxx.in -D BLOCKSIZES="${BLOCKSIZES}" -D outputfile=${CMAKE_CURRENT_BINARY_DIR}/amgcl_dl.cxx  -D TVTI="DL"  -D TV="double"  -D TI="long"  -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/specialize.cmake
    DEPENDS amgcl_tv_ti.cxx.in ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/specialize.cmake  amgcl_c_dl.h
)

add_custom_command(
    OUTPUT  block_amgcl_di.cxx 
    COMMAND ${CMAKE_COMMAND} -D inputfile=${CMAKE_CURRENT_SOURCE_DIR}/block_amgcl_tv_ti.cxx.in -D BLOCKSIZES="${BLOCKSIZES}" -D outputfile=${CMAKE_CURRENT_BINARY_DIR}/block_amgcl_di.cxx  -D TVTI="DI" -D TV="double"  -D TI="int"  -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/specialize.cmake
    DEPENDS block_amgcl_tv_ti.cxx.in ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/specialize.cmake  amgcl_c_di.h
)

add_custom_command(
    OUTPUT  block_amgcl_dl.cxx
    COMMAND ${CMAKE_COMMAND} -D inputfile=${CMAKE_CURRENT_SOURCE_DIR}/block_amgcl_tv_ti.cxx.in -D BLOCKSIZES="${BLOCKSIZES}" -D outputfile=${CMAKE_CURRENT_BINARY_DIR}/block_amgcl_dl.cxx  -D TVTI="DL"  -D TV="double"  -D TI="long"  -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/specialize.cmake
    DEPENDS block_amgcl_tv_ti.cxx.in ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/specialize.cmake  amgcl_c_dl.h
)


add_library(amgcl_c SHARED ${SOURCES} default_params.cxx)
add_dependencies(amgcl_c amgcl_c.h)
target_include_directories(amgcl_c
                           PUBLIC
                           "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/../include"
                           "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>/../include"
                           "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_include_directories(amgcl_c SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})
        

