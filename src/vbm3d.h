/*
 * Copyright (c) 2018, Thibaud Ehret <ehret.thibaud@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef VBM3D_H_INCLUDED
#define VBM3D_H_INCLUDED

#include <fftw3.h>
#include <vector>
#include "Utilities/LibVideoT.hpp"
#include "parameters.h"

//#define OPTICALFLOW

/** ------------------ **/
/** - Main functions - **/
/** ------------------ **/
//! Main function
int run_vbm3d(
    const float sigma
,   Video<float> &img_noisy
#ifdef OPTICALFLOW
,   Video<float> &flow
#endif
,   Video<float> &img_basic
,   Video<float> &img_denoised
,   const Parameters& prms_1
,   const Parameters& prms_2
,   const unsigned color_space
);

//! 1st step of VBM3D
void vbm3d_1st_step(
    const float sigma
,   Video<float> const& img_noisy
#ifdef OPTICALFLOW
,   Video<float> &flow
#endif
,   Video<float> &img_basic
,   const Parameters& prms
,   fftwf_plan *  plan_2d
,   fftwf_plan *  plan_2d_inv
,   VideoUtils::CropPosition* crop
,   const unsigned color_space
,   Video<float>& originalVideo
,   Video<float>& numberator
,   Video<float>& denominator
);

//! 2nd step of VBM3D
void vbm3d_2nd_step(
    const float sigma
,   Video<float> const& img_noisy
,   Video<float> const& img_basic
#ifdef OPTICALFLOW
,   Video<float> &flow
#endif
,   Video<float> &img_denoised
,   const Parameters& prms
,   fftwf_plan *  plan_2d
,   fftwf_plan *  plan_2d_inv
,   VideoUtils::CropPosition* crop
,   const unsigned color_space
,   Video<float>& originalVideo_noisy
,   Video<float>& originalVideo_basic
,   Video<float>& numberator
,   Video<float>& denominator
);

//! Process 2D dct of a group of patches
void dct_2d_process(
    std::vector<float> &DCT_table_2D
,   Video<float> const& img
,   std::vector<unsigned> const& patch_table 
,   fftwf_plan * plan
,   const unsigned kHW
,   const unsigned ktHW
,   std::vector<float> const& coef_norm
);

int computeSimilarPatches(
	std::vector<float>& distances
,	std::vector<unsigned>& indexes
,	unsigned idx
,	const Video<float>& vid
#ifdef OPTICALFLOW
,   Video<float> &flow
#endif
,	const Parameters& prms
);

//! Process 2D bior1.5 transform of a group of patches
void bior_2d_process(
    std::vector<float> &bior_table_2D
,   Video<float> const& img
,   std::vector<unsigned> const& patch_table 
,   const unsigned nHW
,   const unsigned kHW
,   const unsigned ktHW
,   std::vector<float> &lpd
,   std::vector<float> &hpd
);

void dct_2d_inv(
    std::vector<float> &group_3D_table
,   const unsigned kHW
,   const unsigned ktHW
,   const unsigned N
,   std::vector<float> const& coef_norm_inv
,   fftwf_plan * plan
);

void bior_2d_inv(
    std::vector<float> &group_3D_table
,   const unsigned kHW
,   const unsigned ktHW
,   std::vector<float> const& lpr
,   std::vector<float> const& hpr
);

//! HT filtering using Welsh-Hadamard transform (do only
//! third dimension transform, Hard Thresholding
//! and inverse Hadamard transform)
void ht_filtering_hadamard(
    std::vector<float> &group_3D
,   std::vector<float> &tmp
,   const unsigned nSx_r
,   const unsigned kHard
,   const unsigned ktHard
,   const unsigned chnls
,   std::vector<float> const& sigma_table
,   const float lambdaThr3D
,   std::vector<float> &weight_table
);

//! HT filtering using Haar transform (do only
//! third dimension transform, Hard Thresholding
//! and inverse Hadamard transform)
void ht_filtering_haar(
    std::vector<float> &group_3D
,   std::vector<float> &tmp
,   const unsigned nSx_r
,   const unsigned kHard
,   const unsigned ktHard
,   const unsigned chnls
,   std::vector<float> const& sigma_table
,   const float lambdaThr3D
,   std::vector<float> &weight_table
);

//! Wiener filtering using Welsh-Hadamard transform
void wiener_filtering_hadamard(
    std::vector<float> &group_3D_img
,   std::vector<float> &group_3D_est
,   std::vector<float> &tmp
,   const unsigned nSx_r
,   const unsigned kWien
,   const unsigned ktWien
,   const unsigned chnls
,   std::vector<float> const& sigma_table
,   std::vector<float> &weight_table
);

//! Wiener filtering using Haar transform
void wiener_filtering_haar(
    std::vector<float> &group_3D_img
,   std::vector<float> &group_3D_est
,   std::vector<float> &tmp
,   const unsigned nSx_r
,   const unsigned kWien
,   const unsigned ktWien
,   const unsigned chnls
,   std::vector<float> const& sigma_table
,   std::vector<float> &weight_table
);

//! Apply a bior1.5 spline wavelet on a vector of size N x N.
void bior1_5_transform(
    std::vector<float> const& input
,   std::vector<float> &output
,   const unsigned N
,   std::vector<float> const& bior_table
,   const unsigned d_i
,   const unsigned d_o
,   const unsigned N_i
,   const unsigned N_o
);

void temporal_transform(
    std::vector<float>& group_3D
,   const unsigned kHW
,   const unsigned ktHW
,   const unsigned chnls
,   const unsigned nSx_r
);

void temporal_inv_transform(
std::vector<float>& group_3D
,   const unsigned kHW
,   const unsigned ktHW
,   const unsigned chnls
,   const unsigned nSx_r
);

/** ---------------------------------- **/
/** - Preprocessing / Postprocessing - **/
/** ---------------------------------- **/
//! Preprocess coefficients of the Kaiser window and normalization coef for the DCT
void preProcess(
    std::vector<float> &kaiserWindow
,   std::vector<float> &coef_norm
,   std::vector<float> &coef_norm_inv
,   const unsigned kHW
);

#endif // VBM3D_H_INCLUDED
