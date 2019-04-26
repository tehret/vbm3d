#!/bin/bash
#set -x
# dctdenoising_multi.sh 
#        noisy.tif
#        40
#        out.tif      # also creates out.tif.non.tif
#        "-w 8 -1"    # dctdenoising params (optional)
#        3            # LEVELS of pyramid (default: -1 = auto) (optional)
#        2            # R_PYR pyramid ratio: 2 (optional)
#        0.7          # PAR_PYR recomposition ratio : 0.7 (optional)

PATH_MULTISCALE=./

if [ $# -lt 4 ]; then
    echo "$0 noisy.tif sigma out.tif"
    echo "    \"-w 8 -1\"  # dctdenoising params (optional)"
    echo "    -1           # LEVELS of pyramid (default: -1 = auto) (optional)"
    echo "    2            # R_PYR pyramid ratio: 2 (optional), 1.5 also possible"
    echo "    0.7          # PAR_PYR recomposition ratio : 0.7 (optional)"
    echo "    1 150        # first_frame and last_frame for the video"
    echo "    1            # 3D dct for the video instead of frame by frame"
    exit 1
fi

INPUT=$1
NOISE=$2
OUT=$3
OUTPUT=$4
DEN_ARGS=$5
LEVELS=-1
FIRST=1
LAST=150
R_PYR=2   
PAR_PYR=0.7


if [ -n "$6" ]; then
    LEVELS=$6
fi
if [ -n "$7" ]; then
    FIRST=$7
fi
if [ -n "$8" ]; then
    LAST=$8
fi
if [ -n "$9" ]; then
    R_PYR=$9
fi
if [ -n "${10}" ]; then
    PAR_PYR=${10}
fi

# determine the levels based on the image size
if [ $LEVELS -eq -1 ]; then
    PIXELS=$(./num_pixels $INPUT)
    echo $PIXELS
    if [ ${PIXELS} -lt 500000 ]; then
        LEVELS=1
    elif [ ${PIXELS} -lt 2000000 ]; then
        LEVELS=2
    elif [ ${PIXELS} -lt 8000000 ]; then
        LEVELS=3
    else
        LEVELS=4
    fi
    echo "Scales: $LEVELS (auto)"
else
    echo "Scales: $LEVELS (requested)"
fi


DEC_ARGS="-r ${R_PYR}"
MS_ARGS="-c ${PAR_PYR}"

# Clean noisy
rm noisy/*
rm levels/*

$PATH_MULTISCALE/addnoise -i ${INPUT} -sigma $NOISE -f ${FIRST} -l ${LAST} -s $LEVELS # Produces noisy_%04d.tiff
#$PATH_MULTISCALE/addnoise -i ${INPUT} -sigma 0 -f ${FIRST} -l ${LAST} -s $LEVELS # Produces noisy_%04d.tiff

for ((frame=FIRST; frame<=LAST; ++frame))
do
    $PATH_MULTISCALE/decompose $(printf "noisy/noisy_%04d.tiff" $frame) levels/level_ ${LEVELS} _$(printf "%04d" $frame).tiff
done

for ((lvl=LEVELS-1; lvl>=0; --lvl))
do
    sigma=$(bc <<< "scale=2; $NOISE / ${R_PYR}^$lvl")
    python denoise_vbm3d.py levels/ denoised/ level_${lvl}_\*.tiff --sigma $sigma $DEN_ARGS
done

wait
for ((frame=FIRST; frame<=LAST; ++frame))
do
    $PATH_MULTISCALE/recompose denoised/level_ ${LEVELS} _$(printf "%04d" $frame).tiff $(printf $OUT/$OUTPUT $frame)
done

# Copy noisy to ouput folder
mv noisy/* $OUT
# Copy first scale to output folder
mv denoised/level_0_* $OUT
