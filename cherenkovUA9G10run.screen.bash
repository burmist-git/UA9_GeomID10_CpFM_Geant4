#!/bin/bash

runID=$1

screenname="g4CpFMG10_"$runID

exeName=$2
macFileName=$3
seed=$4
outputRootFileName=$5
nameOfTheParticle=$6

screen -S $screenname -L -d -m ./cherenkovUA9G10run.bash $exeName $macFileName $seed $outputRootFileName $nameOfTheParticle
