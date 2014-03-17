#!/bin/bash -x

# Run sequential insertion sort with and without LCP

set -e
source tb-20140109-hostsizes.sh

HOST=`hostname -s`
TESTNAME=inssort

mkdir -p output

source ./fresh-build.sh > /dev/null

test_run() {
    local SHORT=$1
    shift

    OUT="output/$TESTNAME-$HOST-$SHORT.txt"

    rm -vf "$OUT"

    OPTIONS="$ALGO -r 3 --no-check --fork -M mmap_node0"

    # expand HOSTSIZE
    HOSTSIZE=`hostsize $SHORT`

    HOSTSIZE=$(( 1024*1024 ))
    SIZE=16
    REPSIZE=$(( 128*1024*1024 ))

    while [[ $SIZE -le $HOSTSIZE ]]; do

        REPEAT=$(( $REPSIZE / $SIZE ))

        # run once with full threads
        $PSSTEST $OPTIONS -s $SIZE -R $REPEAT $@ | tee -a $OUT

        SIZE=$(( $SIZE * 2 ))

    done 
}

ALGO="-A insertion_sort -A bingmann/lcp_insertion_sort -A bs/mkqsort"

test_run randomASCII ~/global_data/strings/randomASCII.34359738368.gz

test_run urls ~/global_data/strings/urls.txt.75914446213.gz

test_run gov2 ~/global_data/strings/gov2.457165206582.gz

test_run enwiki --suffix ~/global_data/strings/enwiki.83339804349.gz

test_run set5_url ~/global_data/strings/sinha/set5_url.dat

test_run set6_genome ~/global_data/strings/sinha/set6_genome.dat

test_run set6_nodup ~/global_data/strings/sinha/set6_nodup.dat
