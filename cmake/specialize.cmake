message(STATUS "specialize")
string(REPLACE "\\" "" myblocksizes  "${BLOCKSIZES}")

message(STATUS ${inputfile})

file(READ ${inputfile} string1)
string(REPLACE "TvTi" "${TVTI}" string2  "${string1}")
string(REPLACE "Tv" "${TV}" string3  "${string2}")
string(REPLACE "Ti" "${TI}" string4  "${string3}")
string(REPLACE "BLOCKSIZES" "${myblocksizes}" string5  "${string4}")
file(WRITE ${outputfile} "${string5}")
