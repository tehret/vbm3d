#!/bin/bash
# Computes tvl1 optical flow for a (noisy) sequence. 

I=${1:-""}
F=${2:-1}
L=${3:-1}
O=${4:-""}

for i in `seq $F $(($L - 1))`;
do
    ./tvl1flow `printf $D $((i + 1))` \
        `printf $D $i` \
        `printf $O $((i + 1))` \
        4 0.25 0.2 0.3 100 2 0.5 5 0.01 0;
done
