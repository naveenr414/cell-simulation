#!/bin/bash

#very basic script to run x simulations with a particular parameter setting

rundir=$1
nrruns=$2
parfile=$3
nrcores=$4

exewrap="~/Cell_Evolution_stickymoves/scripts/exewrap.sh"
exe="./cell_evolution"

if [ -d $rundir ]; then
  echo "Directory $rundir already exists. Exiting..."
  exit 1

else
  mkdir $rundir
  cp $parfile $rundir
  parfile=$rundir/$(basename $parfile)
fi


#turn this into a gnu parallel?
parallel --gnu -j $nrcores $exewrap $exe $parfile $rundir <<< "$(seq $nrruns)"
