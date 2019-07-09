#!/bin/bash
# Script that computes both optical flows and uses them for VBM3D denoising

INPUT=$1
NOISE=$2
OUTPUT=$3
FIRST=$4
LAST=$5
OPTIONS=${6:-""}

./tvl1flow_backward.sh $INPUT $FIRST $LAST backward_%04d.flo
./tvl1flow_forward.sh $INPUT $FIRST $LAST forward_%04d.flo

./VBM3Ddenoising -i $INPUT -deno $OUTPUT -f $FIRST -l $LAST -sigma $NOISE -bflow backward_%04d.flo -fflow forward_%04d.flo $OPTIONS
