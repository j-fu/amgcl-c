message(STATUS "specialize")

message(STATUS ${inputfile})

file(READ ${inputfile} string1)
string(REPLACE "TvTi" "${TVTI}" string2  "${string1}")
string(REPLACE "Tv" "${TV}" string3  "${string2}")
string(REPLACE "Ti" "${TI}" string4  "${string3}")
file(WRITE ${outputfile} "${string4}")
