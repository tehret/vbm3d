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

# compute backward optical flow {{{1
BFLO="$OUT/bflo-%04d.flo"
FFLO="$OUT/fflo-%04d.flo"

# configure optical flow method
read -ra O <<< "$OPM"
OFBIN="$DIR/${O[0]}"

NPROC=2
case ${O[0]} in
	"tvl1flow")
		FSCALE=${O[1]}; DW=${O[2]};
		OFPRMS="$NPROC 0 $DW 0 0 $FSCALE";;
		# nproc tau lambda theta nscales fscale zfactor nwarps epsilon verbos
	"phsflow")
		FSCALE=${O[1]}; ALPHA=${O[2]};
		OFPRMS="$NPROC $ALPHA 0 $FSCALE";;
		# nproc alpha nscales fscale zfactor nwarps TOL maxiter verbose
	"rof")
		FSCALE=${O[1]}; ALPHA=${O[2]}; GAMMA=${O[3]};
		OFPRMS="$NPROC $ALPHA $GAMMA 10 $FSCALE 0.5 1e-4 1 8";;
		# nproc alpha gamma nscales fscale zfactor TOL inner_iter outer_iter verbose
	"rdpof")
		FSCALE=${O[1]}; ALPHA=${O[2]}; GAMMA=${O[3]};
		OFPRMS="$NPROC 3 $ALPHA $GAMMA 0 10 $FSCALE 0.5 1e-4 1 8";;
		# nproc method alpha gamma lambda nscales fscale zfactor TOL i_iter o_iter verbose
	* )
		echo ERROR: unknown optical flow $OFBIN
		exit 1;;
esac

# backward flow
for i in $(seq $((FFR+1)) $LFR);
do
	file=$(printf $BFLO $i)
	if [ ! -f $file ]; then
		$OFBIN $(printf $SEQ $i) $(printf $SEQ $((i-1))) $file $OFPRMS
	fi
done

# forward optical flow
for i in $(seq $((FFR)) $((LFR-1)));
do
	file=$(printf $FFLO $i)
	if [ ! -f $file ]; then
		$OFBIN $(printf $SEQ $i) $(printf $SEQ $((i+1))) $file $OFPRMS
	fi
done


# denoise sequence {{{1
$DIR/VBM3Ddenoising -i $SEQ -f $FFR -l $LFR -sigma $SIG -add 0 \
	-bflow $OUT/bflo-%04d.flo -fflow $OUT/fflo-%04d.flo $DPM \
	-deno $OUT/deno-%04d.tiff

# vim:set foldmethod=marker:
