#!/bin/bash -x

# Run selection of sequential algorithms for paper

set -e

HOST=`hostname -s`
TESTNAME=seqalgos-memuse

mkdir -p output

source ./fresh-build.sh -DWITH_MALLOC_COUNT=ON > /dev/null

MiB=2**20
GiB=2**30

HOSTSIZE_urls=$(( 4 * $GiB ))
HOSTSIZE_randomASCII=$(( 4 * $GiB ))
HOSTSIZE_gov2=$(( 4 * $GiB ))
HOSTSIZE_enwiki=$(( 256 * $MiB ))
HOSTSIZE_set5_url=$(( 999 * $MiB ))
HOSTSIZE_set6_genome=$(( 999 * $MiB ))
HOSTSIZE_set6_nodup=$(( 999 * $MiB ))

TIMEOUT_urls=900
TIMEOUT_randomASCII=1800
TIMEOUT_gov2=1200
TIMEOUT_enwiki=1200
TIMEOUT_set5_url=900
TIMEOUT_set6_genome=900
TIMEOUT_set6_nodup=900

test_run() {
    local SHORT=$1
    shift

    OUT="output/$TESTNAME-$HOST-$SHORT.txt"
    #rm -f "$OUT"
    [ -e "$OUT" ] && return

    OPTIONS="-r 3 --no-check --datafork "
    ALGOS=""
    ALGOS="$ALGOS -A bingmann/msd_CI5 "
    ALGOS="$ALGOS -A rantala/msd_CE7 "
    ALGOS="$ALGOS -A bingmann/sample_sortBTCU "
    ALGOS="$ALGOS -A bingmann/sample_sortBTCE3U "
    ALGOS="$ALGOS -A bingmann/sequential_mkqs_cache8 "
    ALGOS="$ALGOS -A bingmann/stdsort1 "
    ALGOS="$ALGOS -A bs/mkqsort "
    ALGOS="$ALGOS -A ng/cradix "
    ALGOS="$ALGOS -A ng/lcpmergesort "
    ALGOS="$ALGOS -A sinha/burstsortA "
    ALGOS="$ALGOS -A sinha/fbC_burstsort "
    ALGOS="$ALGOS -A sinha/sCPL_burstsort "

    # expand HOSTSIZE
    local VAR=HOSTSIZE_${SHORT}
    HOSTSIZE=${!VAR:?Error HOSTSIZE is not defined!}

    # expand HOSTSIZE
    VAR=TIMEOUT_${SHORT}
    TIMEOUT=${!VAR:?Error TIMEOUT is not defined!}

    # run once all basic algorithms with interleaved string memory
    $PSSTEST $OPTIONS $ALGOS -T $TIMEOUT -s $HOSTSIZE $@ | tee -a "$OUT"
}

test_run set5_url ~/global_data/strings/sinha/set5_url.dat

test_run set6_genome ~/global_data/strings/sinha/set6_genome.dat

test_run set6_nodup ~/global_data/strings/sinha/set6_nodup.dat

test_run gov2 ~/global_data/strings/gov2.457165206582.lzo

test_run urls ~/global_data/strings/urls.txt.75914446213.gz

test_run randomASCII ~/global_data/strings/randomASCII.34359738368.lzo

test_run enwiki --suffix ~/global_data/strings/enwiki.83339804349.gz
