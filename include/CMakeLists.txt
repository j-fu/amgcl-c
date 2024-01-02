add_custom_target(amgcl_c_di.h
    ${CMAKE_COMMAND} -D inputfile=${CMAKE_CURRENT_SOURCE_DIR}/amgcl_c_tv_ti.h.in -D BLOCKSIZES="${BLOCKSIZES}" -D outputfile=${CMAKE_CURRENT_BINARY_DIR}/amgcl_c_di.h  -D TVTI="DI" -D TV_TI="di" -D TV="double"  -D TI="int"  -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/specialize.cmake
    DEPENDS amgcl_c_tv_ti.h.in ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/specialize.cmake
)

add_custom_target(amgcl_c_dl.h
    ${CMAKE_COMMAND} -D inputfile=${CMAKE_CURRENT_SOURCE_DIR}/amgcl_c_tv_ti.h.in -D BLOCKSIZES="${BLOCKSIZES}" -D outputfile=${CMAKE_CURRENT_BINARY_DIR}/amgcl_c_dl.h  -D TVTI="DL"  -D TV_TI="dl" -D TV="double"  -D TI="long"  -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/specialize.cmake
    DEPENDS amgcl_c_tv_ti.h.in ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/specialize.cmake
)

add_custom_target(amgcl_c.h
  ${CMAKE_COMMAND} -E cat ${HEADERS} > ${CMAKE_BINARY_DIR}/include/amgcl_c.h
  DEPENDS ${HEADERS}
  )