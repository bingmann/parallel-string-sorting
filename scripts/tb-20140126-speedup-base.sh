#!/bin/bash -x

# Run base sequential tests for speedups for selected algorithms

set -e
source tb-20140109-hostsizes.sh

HOST=`hostname -s`
TESTNAME=speedup-base

mkdir -p output

source ./fresh-build.sh > /dev/null

test_run() {
    local SHORT=$1
    shift

    OUT="output/$TESTNAME-$HOST-$SHORT.txt"
    #rm -f "$OUT"
    [ -e "$OUT" ] && return

    OPTIONS="-r 2 --no-check --fork -M mmap "
    ALGOS=""
    ALGOS="$ALGOS -A bingmann/sequential_mkqs_cache8"
    ALGOS="$ALGOS -A rantala/msd_CE5 "
    ALGOS="$ALGOS -A rantala/msd_CE6 "
    ALGOS="$ALGOS -A rantala/msd_CE7 "
    ALGOS="$ALGOS -A rantala/multikey_cache8 "
    ALGOS="$ALGOS -A rantala/msd_CI_adaptive "

    # expand HOSTSIZE
    hostsize $SHORT

    # run all basic sequential algorithms with plain string memory
    $PSSTEST $OPTIONS $ALGOS -s $HOSTSIZE $@ | tee -a "$OUT"
}

test_run randomASCII ~/global_data/strings/randomASCII.34359738368.gz

test_run urls ~/global_data/strings/urls.txt.75914446213.gz

test_run gov2 ~/global_data/strings/gov2.457165206582.gz

test_run enwiki --suffix ~/global_data/strings/enwiki.83339804349.gz

test_run set5_url ~/global_data/strings/sinha/set5_url.dat

test_run set6_genome ~/global_data/strings/sinha/set6_genome.dat

test_run set6_nodup ~/global_data/strings/sinha/set6_nodup.dat
