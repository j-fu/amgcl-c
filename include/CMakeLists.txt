add_custom_command(
    OUTPUT amgcl_c_di.h
    COMMAND ${CMAKE_COMMAND} -D inputfile=${CMAKE_CURRENT_SOURCE_DIR}/amgcl_c_tv_ti.h.in
                             -D BLOCKSIZES="${BLOCKSIZES}"
                             -D outputfile=${CMAKE_BINARY_DIR}/include/amgcl_c_di.h
                             -D TVTI="DI"
                             -D TV_TI="di"
                             -D TV="double"
                             -D TI="int"
                             -P ${CMAKE_SOURCE_DIR}/cmake/specialize.cmake
                             DEPENDS amgcl_c_tv_ti.h.in
                             )

add_custom_command(
  OUTPUT amgcl_c_dl.h
  COMMAND ${CMAKE_COMMAND} -D inputfile=${CMAKE_CURRENT_SOURCE_DIR}/amgcl_c_tv_ti.h.in
                           -D BLOCKSIZES="${BLOCKSIZES}"
                           -D outputfile=${CMAKE_BINARY_DIR}/include/amgcl_c_dl.h
                           -D TVTI="DL"
                           -D TV_TI="dl"
                           -D TV="double"
                           -D TI="${INT64_T}"
                           -P ${CMAKE_SOURCE_DIR}/cmake/specialize.cmake
                           DEPENDS amgcl_c_tv_ti.h.in
                           )
                         
add_custom_command(
  OUTPUT amgcl_c.h
  COMMAND ${CMAKE_COMMAND} -E cat ${HEADERS} > amgcl_c.h
  DEPENDS ${HEADERS}
  )

add_custom_target(amgcl_c_h DEPENDS amgcl_c.h)
