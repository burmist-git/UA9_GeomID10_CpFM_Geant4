#!/bin/bash

make clean; make -j8;

./cherenkovUA9G10 vis.mac 4313 ../dataSim/cherenkovUA9G10.root proton
