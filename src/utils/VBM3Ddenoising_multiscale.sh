#!/bin/bash
# Script that perform a multiscale denoising using VBM3D

if [ $# -lt 4 ]; then
    echo "$0 input_%03d.png sigma outputFolder out_%03d.tif"
    echo "    1     # Index of the first frame of the video, (optional) default is 1"
    echo "    100   # Index of the last frame of the video, (optional) default is 100"
    echo "    \"\"  # Denoising params. Used to provide arguments to VBM3D similarly to how they would be provided in the command line. For example if one wants to use a patch size of 16 for the first step, this parameter becomes \"-kHard 16\". \"\" keep the default parameters"
    echo "    3     # Number of scales, (optional) default is 3"
    echo "    0     # Type of multiscaler (0: DCT, 1: Gaussian, 2: Lanczos), default is 0 (optional)"
    echo "    0.7   # Recomposition ratio : 0.7 (optional)"
    exit 1
fi

INPUT=$1
NOISE=$2
OUT=$3
OUTPUT=$4
FIRST=${5:-1}
LAST=${6:-100}
DEN_ARGS=${7:-""}
LEVELS=${8:-3}
MSTYPE=${9:-0}
PAR_PYR=${10:-0.7}
R_PYR=2

# We assume that the binaries are in the same folder as the script
DIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

MS_ARGS="-c ${PAR_PYR}"

# Create folders
mkdir -p $OUT/noisy
mkdir -p $OUT/levels
mkdir -p $OUT/denoised

# Add noise
$DIR/addnoise -i ${INPUT} -sigma $NOISE -f ${FIRST} -l ${LAST} \
	-o "$OUT/noisy/noisy_%04d.tiff"

# Decompose the sequence into the multiple scales
for ((frame=FIRST; frame<=LAST; ++frame))
do
    $DIR/decompose $(printf "$OUT/noisy/noisy_%04d.tiff" $frame) \
		 "$OUT/levels/level_" ${LEVELS} "_"$(printf "%04d" $frame)".tiff" -t $MSTYPE
done

# Denoise the different scales
for ((lvl=LEVELS-1; lvl>=0; --lvl))
do
    sigma=$(bc <<< "scale=2; $NOISE / ${R_PYR}^$lvl")
    $DIR/VBM3Ddenoising -i "$OUT/levels/level_${lvl}_%04d.tiff" \
		 -deno "$OUT/denoised/level_${lvl}_%04d.tiff" \
		 -f ${FIRST} -l ${LAST} -add false -sigma $sigma $DEN_ARGS
done

# Recompose the different scales into the final denoised sequence
for ((frame=FIRST; frame<=LAST; ++frame))
do
    $DIR/recompose $OUT/denoised/level_ ${LEVELS} \
		 _$(printf "%04d" $frame).tiff $(printf $OUT/$OUTPUT $frame) \
		 -t $MSTYPE ${MS_ARGS}
done

# Compute the final PSNRs
./psnr -i $OUT/denoised/level_0_%04d.tiff -r ${INPUT} -f ${FIRST} -l ${LAST} > $OUT/psnr-ss
./psnr -i $OUT/$OUTPUT -r ${INPUT} -f ${FIRST} -l ${LAST} > $OUT/psnr-ms
