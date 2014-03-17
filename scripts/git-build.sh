#!/bin/bash -x

set -e

BUILD=~/`hostname -s`/pss-opt

# remove and mkdir $BUILD directory
rm -rvf $BUILD
mkdir -p $BUILD

pushd $BUILD

# fresh cmake configure and make
cmake ~/parallel-string-sorting -DCMAKE_BUILD_TYPE=Release "$@"
make -j || exit

# go into src and run psstest
cd $BUILD/src
pwd
./psstest

export PSSTEST=`pwd`/psstest

popd
