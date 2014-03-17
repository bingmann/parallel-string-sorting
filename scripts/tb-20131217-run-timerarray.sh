#!/bin/bash -x

set -e

HOST=`hostname -s`
TESTNAME=timerarray

# Run TimerArray tests with quadraticly increasing input sizes

source ./fresh-build.sh -DCMAKE_CXX_FLAGS="-DTIMERARRAY_REAL"
echo "Fresh build finished"

mkdir -p output

test_run() {
    local SIZE_LOW=$1
    local SIZE_HIGH=$2
    local INSHORT=$3
    local INPUT=$4
    local EXTRA=$5

    local SIZE=$SIZE_LOW

    while [ $SIZE -le $SIZE_HIGH ]; do

        OUT="output/$TESTNAME-$HOST-$INSHORT.txt"

        OPTIONS="-a bingmann/parallel_sample_sortBTCU2_nolcp -r 3 $EXTRA"

        $PSSTEST $OPTIONS -s "${SIZE}mb" $INPUT | tee -a $OUT

        # double size
        SIZE=$(( 2 * $SIZE ))
    done
}

test_run 16 40000 randomASCII ~/global_data/strings/randomASCII.34359738368.gz

test_run 16 80000 urls ~/global_data/strings/urls.txt.75914446213.gz

test_run 16 500000 gov2 ~/global_data/strings/gov2.457165206582.gz

test_run 16 8000 enwiki ~/global_data/sac-corpus/bigworld/enwiki-20120601-pages-meta-current.xml.83339804349.gz --suffix
