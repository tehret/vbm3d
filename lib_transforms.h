/*
 * Copyright (c) 2011, Marc Lebrun <marc.lebrun@cmla.ens-cachan.fr>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef LIB_TRANSFORMS_INCLUDED
#define LIB_TRANSFORMS_INCLUDED

#include<vector>
#include"Utilities/LibVideoT.hpp"

#define SQRT2     1.414213562373095
#define SQRT2_INV 0.7071067811865475

//! Compute a Bior1.5 2D
void bior_2d_forward(
    Video<float> const& input
,   std::vector<float> &output
,   const unsigned N
,   const unsigned x
,   const unsigned y
,   const unsigned t
,   const unsigned c
,   const unsigned d_o
,   std::vector<float> const& lpd
,   std::vector<float> const& hpd
);

void bior_2d_forward_test(
    std::vector<float> const& input
,   std::vector<float> &output
,   const unsigned N
,   const unsigned d_i
,   const unsigned r_i
,   const unsigned d_o
,   std::vector<float> const& lpd
,   std::vector<float> const& hpd
,   std::vector<float> &tmp
,   std::vector<unsigned> &ind_per
);

//! Compute a Bior1.5 2D inverse
void bior_2d_inverse(
    std::vector<float> &signal
,   const unsigned N
,   const unsigned d_s
,   std::vector<float> const& lpr
,   std::vector<float> const& hpr
);

//! Precompute the Bior1.5 coefficients
void bior15_coef(
    std::vector<float> &lp1
,   std::vector<float> &hp1
,   std::vector<float> &lp2
,   std::vector<float> &hp2
);

//! Apply Walsh-Hadamard transform (non normalized) on a vector of size N = 2^n
void hadamard_transform(
    std::vector<float> &vec
,   std::vector<float> &tmp
,   const unsigned N
,   const unsigned d
);

void haar_forward(
    std::vector<float> &vec
,   std::vector<float> &tmp
,   const unsigned N
,   const unsigned d
);

void haar_inverse(
    std::vector<float> &vec
,   std::vector<float> &tmp
,   const unsigned N
,   const unsigned d
);

//! Process the log2 of N
unsigned log2(
    const unsigned N
);

//! Obtain index for periodic extension
void per_ext_ind(
    std::vector<unsigned> &ind_per
,   const unsigned N
,   const unsigned L
);

#endif // LIB_TRANSFORMS_INCLUDED
