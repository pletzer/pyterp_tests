#!/bin/bash -e

## Usage: collect.sh <NAME>

# module loads
ml Python-Iris VTune

## program to run 
prog="python sigrid_conserve.py"
## arguments to pass to the program
args=
## command to run gprof2dot.py
#gprof2dot=gprof2dot.py
gprof2dot=/projects/nesi99999/ANTS-regridding/local/bin/gprof2dot.py

## check argument
if (( $# < 1 ))
then
    echo "Usage: $0 <UNIQUE_NAME> <ARGS>"
    exit 1
fi
if [ -d "$1.dir" ]
then
    echo "Enter unique name"
    exit 2
fi
if (( $# > 1 ))
then
    args="$2"
fi

## run the program to collect data
amplxe-cl -collect hotspots -result-dir $1.dir -search-dir=$SIGRID_DIR/build/lib.linux-x86_64-2.7/sigrid/ -source-search-dir=$SIGRID_DIR/build/cpp -- $prog $args

## generate gprof-like report
amplxe-cl -report gprof-cc -result-dir $1.dir -format text -report-output $1.txt

## create the call graph
$gprof2dot -f axe $1.txt | dot -Tpng -o $1.png
