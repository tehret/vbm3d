#!/bin/bash
PATH_MULTISCALE=./

if [ $# -lt 4 ]; then
    echo "$0 input_%03d.png sigma outputFolder out_%03d.tif"
    echo "    \"\"  # denoising params"
    echo "    3     # Number of scales, (optional) default is 3"
    echo "    1     # index of the first frame of the video, (optional) default is 1"
    echo "    100   # index of the last frame of the video, (optional) default is 100"
    echo "    0     # type of multiscaler (0: DCT, 1: Gaussian, 2: Lanczos), default is 0 (optional)"
    echo "    2     # R_PYR pyramid ratio: 2 (optional), 1.5 also possible"
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
    MSTYPE=$9
fi
if [ -n "${10}" ]; then
    R_PYR=${10}
fi
if [ -n "${11}" ]; then
    PAR_PYR=${11}
fi

DEC_ARGS="-r ${R_PYR}"
MS_ARGS="-c ${PAR_PYR}"

# Clean noisy
rm -r noisy
rm -r levels
rm -r denoised
mkdir noisy
mkdir levels
mkdir denoised

$PATH_MULTISCALE/addnoise -i ${INPUT} -sigma $NOISE -f ${FIRST} -l ${LAST}

for ((frame=FIRST; frame<=LAST; ++frame))
do
    $PATH_MULTISCALE/decompose $(printf "noisy/noisy_%04d.tiff" $frame) levels/level_ ${LEVELS} _$(printf "%04d" $frame).tiff -t $MSTYPE
done

for ((lvl=LEVELS-1; lvl>=0; --lvl))
do
    sigma=$(bc <<< "scale=2; $NOISE / ${R_PYR}^$lvl")
    ./VBM3Ddenoising_OF.sh levels/level_${lvl}_%04d.tiff $sigma denoised/level_${lvl}_%04d.tiff ${FIRST} ${LAST} $DEN_ARGS
done

wait
for ((frame=FIRST; frame<=LAST; ++frame))
do
    $PATH_MULTISCALE/recompose denoised/level_ ${LEVELS} _$(printf "%04d" $frame).tiff $(printf $OUT/$OUTPUT $frame) -t $MSTYPE
done

./psnr -i $OUT/$OUTPUT -r ${INPUT} -f ${FIRST} -l ${LAST} > $OUT/psnr.txt

# Copy noisy to ouput folder
mv noisy $OUT
# Copy noisy to ouput folder
mv levels $OUT
# Copy first scale to output folder
mv denoised $OUT
