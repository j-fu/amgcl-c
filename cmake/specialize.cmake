message(STATUS "specialize")

message(STATUS ${inputfile})

file(READ ${inputfile} string1)
string(REPLACE "TvTi" "${TVTI}" string2  "${string1}")
string(REPLACE "Tv" "${TV}" string3  "${string2}")
string(REPLACE "Ti" "${TI}" string4  "${string3}")
string(REPLACE "DOX" "DO(2) DO(3) DO(4) DO(5) DO(6) DO(7) DO(8) DO(9) DO(10)" string5  "${string4}")
file(WRITE ${outputfile} "${string5}")
