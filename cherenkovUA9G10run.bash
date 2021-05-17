#!/bin/bash

source /home/leonid/root_v5.34.34/root_v5.34.34-install/bin/thisroot.sh

outdataDir='../dataSim/'

if [ $# -eq 0 ]
then
    make clean; make -j8;
    echo " --> ./cherenkovUA9G10 run.mac 4313 $outdataDir/cherenkovUA9G10.root proton"
    ./cherenkovUA9G10 run.mac 4313 $outdataDir/cherenkovUA9G10.root proton
fi

# to run bunch of simulations with screen
if [ $# -eq 5 ]
then
    exeName=$1
    macFileName=$2
    seed=$3
    outputRootFileName=$4
    nameOfTheParticle=$5
    echo $exeName
    echo $macFileName
    echo $seed
    echo $outputRootFileName
    echo $nameOfTheParticle
    outLogFile=$outdataDir/$outputRootFileName'.log'
    ./$exeName $macFileName $seed $outdataDir/$outputRootFileName $nameOfTheParticle | tee $outLogFile
fi
