#!/bin/bash

if [ $# -eq 0 ]
  then
    echo "No arguments supplied; first argument should be parameter file"
    exit 0
fi

cd ../src; rm -rf backup/ data_film2/ data_cellcount.txt; ./cell_evolution ../data/parameters/$1
