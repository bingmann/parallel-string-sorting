#!/bin/bash

# Maximum dataset sizes for different test hosts

set -e

MiB=2**20
GiB=2**30

HOSTSIZE_i10pc127_urls=$(( 999 * $GiB ))
HOSTSIZE_i10pc127_randomASCII=$(( 32 * $GiB ))
HOSTSIZE_i10pc127_gov2=$(( 128 * $GiB ))
HOSTSIZE_i10pc127_enwiki=$(( 4 * $GiB ))

HOSTSIZE_i10pc126_urls=$(( 999 * $GiB ))
HOSTSIZE_i10pc126_randomASCII=$(( 24 * $GiB ))
HOSTSIZE_i10pc126_gov2=$(( 64 * $GiB ))
HOSTSIZE_i10pc126_enwiki=$(( 2 * $GiB ))

HOSTSIZE_i10pc124_urls=$(( 32 * $GiB ))
HOSTSIZE_i10pc124_randomASCII=$(( 12 * $GiB ))
HOSTSIZE_i10pc124_gov2=$(( 24 * $GiB ))
HOSTSIZE_i10pc124_enwiki=$(( 1 * $GiB ))

HOSTSIZE_i10pc113_urls=$(( 4 * $GiB ))
HOSTSIZE_i10pc113_randomASCII=$(( 2 * $GiB ))
HOSTSIZE_i10pc113_gov2=$(( 4 * $GiB ))
HOSTSIZE_i10pc113_enwiki=$(( 256 * $MiB ))

HOSTSIZE_i10pc112_urls=$(( 4 * $GiB ))
HOSTSIZE_i10pc112_randomASCII=$(( 2 * $GiB ))
HOSTSIZE_i10pc112_gov2=$(( 4 * $GiB ))
HOSTSIZE_i10pc112_enwiki=$(( 256 * $MiB ))

HOSTSIZE_i10pc150_urls=$(( 8 * $GiB ))
HOSTSIZE_i10pc150_randomASCII=$(( 3 * $GiB ))
HOSTSIZE_i10pc150_gov2=$(( 8 * $GiB ))
HOSTSIZE_i10pc150_enwiki=$(( 512 * $MiB ))

HOSTSIZE_i10pc151_urls=$(( 8 * $GiB ))
HOSTSIZE_i10pc151_randomASCII=$(( 3 * $GiB ))
HOSTSIZE_i10pc151_gov2=$(( 8 * $GiB ))
HOSTSIZE_i10pc151_enwiki=$(( 512 * $MiB ))

hostsize() {
    local SHORT=$1

    local VAR=HOSTSIZE_${HOST}_${SHORT}

    #HOSTSIZE=100mb
    #return 0

    if [[ -n "${!VAR}" ]]; then
        HOSTSIZE=${!VAR:?Error HOSTSIZE is not defined!}
        return 0
    fi

    if [[ "$SHORT" == "set5_url" || "$SHORT" == "set6_genome" || "$SHORT" == "set6_nodup"  ]]; then
        HOSTSIZE=$(( $GiB ))
        return 0
    fi

    echo "Cannot determine host's $HOST maximum file size for $SHORT"
    exit
}
