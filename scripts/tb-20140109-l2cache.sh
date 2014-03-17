#!/bin/bash -x

# Run first parallel sample sort step with different classification tree sizes.

set -e
source tb-20131222-host-sizes.sh

HOST=`hostname -s`
TESTNAME=l2cache

mkdir -p output

source ./fresh-build.sh -DCMAKE_CXX_FLAGS:STRING="-DPS5_SINGLE_STEP=true" > /dev/null

test_run() {
    local SHORT=$1
    shift

    OUT="output/$TESTNAME-$HOST-$SHORT-$ALGO.txt"

    OPTIONS="-a bingmann/parallel_sample_sort$ALGO -r 3 --no-check --fork -M mmap_node0"

    # expand HOSTSIZE
    HOSTSIZEVAR=HOSTSIZE_${HOST}_${SHORT}
    HOSTSIZE=${!HOSTSIZEVAR:?Error HOSTSIZE is not defined!}

    # run once with full threads
    $PSSTEST $OPTIONS -s $HOSTSIZE $@ --thread-list 16 | tee -a $OUT

    # run again with just one thread
    #$PSSTEST $OPTIONS -s $HOSTSIZE --thread-list 1 $@ | tee -a $OUT
}

for ALGO in BTCU2_nolcp BTCU1_nolcp BTCEU1_nolcp BTCTU2_nolcp BTCE_nolcp BTC_nolcp BSC_nolcp BTCTU1_nolcp BTCT_nolcp; do

    test_run randomASCII ~/global_data/strings/randomASCII.34359738368.gz

    #test_run urls ~/global_data/strings/urls.txt.75914446213.gz

    #test_run gov2 ~/global_data/strings/gov2.457165206582.gz

    #test_run enwiki --suffix ~/global_data/strings/enwiki.83339804349.gz

done
