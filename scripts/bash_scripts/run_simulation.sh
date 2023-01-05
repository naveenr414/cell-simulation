#!/bin/bash

if [ $# -eq 0 ]
  then
    echo "No arguments supplied; first argument should be parameter file, and second the identifier"
    exit 0
fi

backup_folder=backup_$2/
data_folder=data_film_$2/
output_file=data_cellcount_$2.txt
parameter_file=data_cellcount_$2.par
error_file=error_$2.txt

cd ../src; rm -rf $backup_folder $data_folder $output_file
cd ../src; ./cell_evolution ../data/parameters/$1 -datafile $output_file -backupdir $backup_folder -datadir $data_folder > ../src/errors/$error_file
cp ../data/parameters/$1 ../data/output/$parameter_file
cp ../src/$output_file ../data/output/$output_File
