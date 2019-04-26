/*
 * Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef UTILITIES_H_INCLUDED
#define UTILITIES_H_INCLUDED

#include <vector>
#include "LibImages.h"

/**
 * @brief Convenient function to use the sort function provided by the vector library.
 **/
bool comparaisonFirst(
	const std::pair<float, unsigned> &i_pair1
,	const std::pair<float, unsigned> &i_pair2
);

bool comparaisonInverseFirst(
	const std::pair<float, unsigned> &i_pair1
,	const std::pair<float, unsigned> &i_pair2
);

/**
 * @brief Clip a value between min and max
 *
 * @param i_value: value to clip;
 * @param i_min: minimum value;
 * @param i_max: maximum value.
 *
 * @return value clipped between [min, max].
 **/
float clip(
	const float i_value
,	const float i_min
,	const float i_max
);

/**
 * @brief Obtain and substract the baricenter of io_group3d.
 *
 * @param io_group3d(p_rows x p_cols) : data to center;
 * @param o_baricenter(p_cols): will contain the baricenter of io_group3d;
 * @param p_rows, p_cols: size of io_group3d.
 *
 * @return none.
 **/
void centerData(
	std::vector<float> &io_group3d
,	std::vector<float> &o_baricenter
,	const unsigned p_rows
,	const unsigned p_cols
);

/**
 * @brief Compute the average standard deviation of a set of patches.
 *
 * @param i_Set(p_sP, p_nSimP): set of patches;
 * @param p_sP : size of a patch;
 * @param p_nSimP: number of patches in the set;
 * @param p_nChannels: number of channels of the image.
 *
 * @return the average standard deviation of the set
 **/
float computeStdDeviation(
	std::vector<float> const& i_Set
,	const unsigned p_sP
,	const unsigned p_nSimP
,	const unsigned p_nChannels
);

/**
 * @brief Determine a and b such that : n = a * b, with a and b as greatest as possible
 *
 * @param i_n : number to decompose;
 * @param o_a : will contain a;
 * @param o_b : will contain b.
 *
 * @return none.
 **/
void determineFactor(
    const unsigned i_n
,   unsigned &o_a
,   unsigned &o_b
);

/**
 * @brief Write PSNR and RMSE in a .txt for both basic and denoised images.
 *
 * @param p_pathName: path name of the file;
 * @param p_sigma: value of the noise;
 * @param p_psnr: value of the PSNR of the denoised image;
 * @param p_rmse: value of the RMSE of the denoised image;
 * @param p_truncateFile: if true, erase the file when open it. Otherwise
 *        write at the end of the file;
 * @param p_app: in order to specify the image.
 *
 * @return EXIT_FAILURE if the file can't be opened.
 **/
int writingMeasures(
    const char* p_pathName
,   const float p_sigma
,   const float p_psnr
,   const float p_rmse
,   const float p_grps
,   const float p_time
,   const float p_cons
,   const bool  p_truncateFile
,   const char* p_app
);

//! Check if a number is a power of 2
bool power_of_2(
    const unsigned n
);

//! Look for the closest power of 2 number
int closest_power_of_2(
    const unsigned n
);

//! Estimate sigma on each channel according to the choice of the color_space
int estimate_sigma(
    const float sigma
,   std::vector<float> &sigma_table
,   const unsigned chnls
,   const unsigned color_space
);

//! Initialize a 2D fftwf_plan with some parameters
void allocate_plan_2d(
    fftwf_plan* plan
,   const unsigned N
,   const fftwf_r2r_kind kind
,   const unsigned nb
);

//! Initialize a 1D fftwf_plan with some parameters
void allocate_plan_1d(
    fftwf_plan* plan
,   const unsigned N
,   const fftwf_r2r_kind kind
,   const unsigned nb
);

//! Initialize a set of indices
void ind_initialize(
    std::vector<unsigned> &ind_set
,   const unsigned beginning
,   const unsigned end
,   const unsigned step
);

void ind_initialize2(
    std::vector<unsigned> &ind_set
,   const unsigned max_size
,   const unsigned N
,   const unsigned step
);

//! For convenience
unsigned ind_size(
    const unsigned beginning
,   const unsigned end
,   const unsigned step
);


#endif // UTILITIES_H_INCLUDED
