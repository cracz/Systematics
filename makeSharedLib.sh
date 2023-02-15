#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: ./makeSharedLib.sh [file name without .cxx]"
    echo
    exit
fi

FILENAME=$1

if [ ! -e "${FILENAME}.cxx" ]; then
    echo "${FILENAME}.cxx not found here!"
    echo
    exit
else
    genreflex ${FILENAME}.h --rootmap ${FILENAME}.rootmap --rootmap-lib lib${FILENAME}.so
    g++ -o lib${FILENAME}.so ${FILENAME}_rflx.cpp ${FILENAME}.cxx -shared -fPIC `root-config --cflags --libs`
fi
