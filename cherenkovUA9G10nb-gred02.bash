#!/bin/bash

source cleannb-gred02.bash

. setEnvG4WORKDIR.bash

make -j8

#./bin/Linux-g++/cherenkovUA9G10 vis.mac 123123 cherenkovUA9G10.root proton
./bin/Linux-g++/cherenkovUA9G10 run.mac 123123 cherenkovUA9G10.root proton
