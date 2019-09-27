#!/bin/bash
# Runs nlkalman filtering frame by frame

SEQ=$1 # filtered sequence path
FFR=$2 # first frame
LFR=$3 # last frame
SIG=$4 # noise standard dev.
OUT=$5 # output folder
DPM=$6 # denoising parameters
OPM=${7:-"tvl1flow 1 0.40"} # optical flow parameters

# we assume that the binaries are in the same folder as the script
DIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

# error checking {{{1
for i in $(seq $FFR $LFR);
do
	file=$(printf $SEQ $i)
	if [ ! -f $file ]
	then
		echo ERROR: $file not found
		exit 1
	fi
done

mkdir -p $OUT
$DIR/addnoise -i $SEQ -sigma $SIG -f $FFR -l $LFR -o "$OUT/%04d.tiff"
$DIR/vbm3d.sh "$OUT/%04d.tiff" $FFR $LFR $SIG $OUT "$DPM" "$OPM"
$DIR/psnr -i "$OUT/deno-%04d.tiff" -r $SEQ -f $FFR -l $LFR > $OUT/psnr.txt

# vim:set foldmethod=marker:
