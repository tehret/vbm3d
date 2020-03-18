#!/bin/bash

if [ $# -lt 4 ]; then
    echo "$0 input_%03d.png sigma outputFolder out_%03d.tif"
    echo "    \"\"  # denoising params"
    echo "    3     # Number of scales, (optional) default is 3"
    echo "    1     # index of the first frame of the video, (optional) default is 1"
    echo "    100   # index of the last frame of the video, (optional) default is 100"
    echo "    0     # type of multiscaler (0: DCT, 1: Gaussian, 2: Lanczos), default is 0 (optional)"
    echo "    0.7   # PAR_PYR recomposition ratio : 0.7 (optional)"
    exit 1
fi

INPUT=$1
NOISE=$2
OUT=$3
OUTPUT=$4
DEN_ARGS=$5
LEVELS=3
FIRST=1
LAST=100
MSTYPE=0
PAR_PYR=0.7
R_PYR=2

if [ -n "${6}" ];  then LEVELS=${6}  ; fi
if [ -n "${7}" ];  then FIRST=${7}   ; fi
if [ -n "${8}" ];  then LAST=${8}    ; fi
if [ -n "${9}" ];  then MSTYPE=${9}  ; fi
if [ -n "${10}" ]; then PAR_PYR=${10}; fi

# we assume that the binaries are in the same folder as the script
DIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

MS_ARGS="-c ${PAR_PYR}"

# create folders
mkdir -p $OUT/noisy
mkdir -p $OUT/levels
mkdir -p $OUT/denoised

# add noise
$DIR/addnoise -i ${INPUT} -sigma $NOISE -f ${FIRST} -l ${LAST} \
	-o "$OUT/noisy/noisy_%04d.tiff"

for ((frame=FIRST; frame<=LAST; ++frame))
do
    $DIR/decompose $(printf "$OUT/noisy/noisy_%04d.tiff" $frame) \
		 "$OUT/levels/level_" ${LEVELS} "_"$(printf "%04d" $frame)".tiff" -t $MSTYPE
done

for ((lvl=LEVELS-1; lvl>=0; --lvl))
do
    sigma=$(bc <<< "scale=2; $NOISE / ${R_PYR}^$lvl")
    $DIR/VBM3Ddenoising -i "$OUT/levels/level_${lvl}_%04d.tiff" \
		 -deno "$OUT/denoised/level_${lvl}_%04d.tiff" \
		 -f ${FIRST} -l ${LAST} -add false -sigma $sigma $DEN_ARGS
done

for ((frame=FIRST; frame<=LAST; ++frame))
do
    $DIR/recompose $OUT/denoised/level_ ${LEVELS} \
		 _$(printf "%04d" $frame).tiff $(printf $OUT/$OUTPUT $frame) \
		 -t $MSTYPE ${MS_ARGS}
done

./psnr -i $OUT/denoised/level_0_%04d.tiff -r ${INPUT} -f ${FIRST} -l ${LAST} > $OUT/psnr-ss
./psnr -i $OUT/$OUTPUT -r ${INPUT} -f ${FIRST} -l ${LAST} > $OUT/psnr-ms

# convert to png and remove tiffs
 
## # Copy noisy to ouput folder
## mv noisy $OUT
## # Copy noisy to ouput folder
## mv levels $OUT
## # Copy first scale to output folder
## mv denoised $OUT
