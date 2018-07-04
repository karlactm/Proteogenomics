#!/bin/bash
BASE=$1

PATHPEPTRIP="${BASE}/output_peptide_trip/"
if [ ! -d "${PATHPEPTRIP}" ]; then
    mkdir ${PATHPEPTRIP}
fi

echo ARG = ${BASE}
./peptide_trip.pl ${BASE}
