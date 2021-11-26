#!/bin/bash
# Script that computes both optical flows and uses them for VBM3D denoising

if [ $# -lt 3 ]; then
    echo "$0 input_%03d.png sigma out_%03d.tif"
    echo "    1     # Index of the first frame of the video, (optional) default is 1"
    echo "    100   # Index of the last frame of the video, (optional) default is 100"
    echo "    \"\"  # Denoising params. Used to provide arguments to VBM3D similarly to how they would be provided in the command line. For example if one wants to use a patch size of 16 for the first step, this parameter becomes \"-kHard 16\". \"\" keep the default parameters"
    exit 1
fi
INPUT=$1
NOISE=$2
OUTPUT=$3
FIRST=${4:-1}
LAST=${5:-100}
OPTIONS=${6:-""}

# Compute forward and backward flows
./tvl1flow_backward.sh $INPUT $FIRST $LAST backward_%04d.flo
./tvl1flow_forward.sh $INPUT $FIRST $LAST forward_%04d.flo

# Denoise the sequence using the previous two optical flows
./VBM3Ddenoising -i $INPUT -deno $OUTPUT -f $FIRST -l $LAST -sigma $NOISE -bflow backward_%04d.flo -fflow forward_%04d.flo $OPTIONS
