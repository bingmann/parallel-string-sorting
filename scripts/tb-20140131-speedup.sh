#!/bin/bash -x

# Run speedup tests for selected algorithms

set -e
source tb-20140109-hostsizes.sh

HOST=`hostname -s`
TESTNAME=speedup

mkdir -p output

source ./fresh-build.sh > /dev/null

test_run1() {
    local SHORT=$1
    shift

    OUT="output/$TESTNAME-$HOST-$SHORT.txt"
    #rm -f "$OUT"
    [ -e "$OUT" ] && return

    OPTIONS="-r 1 --no-check --fork --some-threads -M mmap_interleave "
    ALGOS=""
    ALGOS="$ALGOS -A bingmann/parallel_sample_sortBTCU2 "
    ALGOS="$ALGOS -A bingmann/parallel_sample_sortBTCEU1 "
    ALGOS="$ALGOS -A bingmann/parallel_mkqs "
    ALGOS="$ALGOS -A bingmann/parallel_radix_sort_8bit "
    ALGOS="$ALGOS -A bingmann/parallel_radix_sort_16bit "
    ALGOS="$ALGOS -A akiba/parallel_radix_sort "

    if [[ "$HOST" == "i10pc112" || "$HOST" == "i10pc151" ]]; then
        ALGOS="$ALGOS -A shamsundar/lcp-merge-string-sort "
        ALGOS="$ALGOS -A rantala/multikey_simd_parallel4 "
        ALGOS="$ALGOS -A rantala/mergesort_lcp_2way_parallel "
        OPTIONS="$OPTIONS -r 15"
    fi

    # expand HOSTSIZE
    hostsize $SHORT

    # run once all basic algorithms with interleaved string memory
    $PSSTEST $OPTIONS $ALGOS -s $HOSTSIZE $@ | tee -a "$OUT"
}

test_run2() {
    local SHORT=$1
    shift

    OUT="output/$TESTNAME-$HOST-$SHORT-segment.txt"
    #rm -f "$OUT"
    [ -e "$OUT" ] && return

    OPTIONS="-r 3 --no-check --fork --some-threads -M mmap_segment "
    ALGOS=""
    ALGOS="$ALGOS -A eberle/ps5-parallel-toplevel-merge "
    #ALGOS="$ALGOS -A eberle/ps5-parallel-toplevel-merge-assisting "

    # expand HOSTSIZE
    hostsize $SHORT

    # run NUMA algorithms with segmented string memory
    $PSSTEST $OPTIONS $ALGOS -s $HOSTSIZE $@ | tee -a "$OUT"
}

test_run() {
    test_run1 "$@"
    test_run2 "$@"
}

test_run set5_url ~/global_data/strings/sinha/set5_url.dat

test_run set6_genome ~/global_data/strings/sinha/set6_genome.dat

test_run set6_nodup ~/global_data/strings/sinha/set6_nodup.dat

test_run gov2 ~/global_data/strings/gov2.457165206582.lzo

test_run urls ~/global_data/strings/urls.txt.75914446213.gz

test_run randomASCII ~/global_data/strings/randomASCII.34359738368.lzo

test_run enwiki --suffix ~/global_data/strings/enwiki.83339804349.gz
