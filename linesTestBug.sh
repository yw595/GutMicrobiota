#!/bin/bash

head -n 2278 "$1/Copcoccus2.faa" | tail -n 18 > "$1/test1.faa"
head -n 2280 "$1/Copcoccus2.faa" | tail -n 20 > "$1/test2.faa"
head -n 2278 "$1/Copcoccus2.faa" | tail -n 20 > "$1/test3.faa"
head -n 2280 "$1/Copcoccus2.faa" | tail -n 22 > "$1/test4.faa"

# all others work, only test4 fails
# note that alone, it is header of 268, and works, must add prev 10 sequences to get error??
# however, adding one letter to get 269 also results error
makeblastdb -in "$1/test4.faa" -dbtype prot -parse_seqids -out "$1/test" -title test
