message(STATUS "specialize ${inputfile} --> ${outputfile}")

string(REPLACE "\\" "" myblocksizes  "${BLOCKSIZES}")

file(READ ${inputfile} string1)
string(REPLACE "TvTi" "${TVTI}" string2  "${string1}")
string(REPLACE "Tv" "${TV}" string3  "${string2}")
string(REPLACE "Ti" "${TI}" string4  "${string3}")
string(REPLACE "BLOCKSIZES" "${myblocksizes}" string5  "${string4}")
file(WRITE ${outputfile} "${string5}")
