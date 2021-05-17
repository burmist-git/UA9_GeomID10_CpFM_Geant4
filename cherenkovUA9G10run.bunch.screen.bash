#!/bin/bash

make clean; make -j8;

nn1=8
nn2=$(( $nn1-1 ))

inexeName='cherenkovUA9G10'
macFileName='run.mac'
nameOfTheParticle='proton'
for i in `seq 0 $nn2`;
do
    runID=$i
    exeName=$inexeName'_'$i
    #echo $exeName
    cp $inexeName $exeName
    seed=$RANDOM
    outputRootFileName=$inexeName'-'$i'.root'
    #echo $runID $exeName $macFileName $seed $outputRootFileName $nameOfTheParticle
    ./cherenkovUA9G10run.screen.bash $runID $exeName $macFileName $seed $outputRootFileName $nameOfTheParticle
done
