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

/**
 * @file utilities.cpp
 * @brief Utilities functions.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

#include "Utilities.h"

#include <math.h>
#include <omp.h>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>

#define YUV       0
#define YCBCR     1
#define OPP       2
#define RGB       3

using namespace std;

/**
 * @brief Convenient function to use the sort function provided by the vector library.
 **/
bool comparaisonFirst(
	const pair<float, unsigned> &i_pair1
,	const pair<float, unsigned> &i_pair2
){
	return i_pair1.first < i_pair2.first;
}

/**
 * @brief Convenient function to use the sort function provided by the vector library. Used to ocmpute the inverse order
 **/
bool comparaisonInverseFirst(
	const pair<float, unsigned> &i_pair1
,	const pair<float, unsigned> &i_pair2
){
	return i_pair1.first > i_pair2.first;
}

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
){
	return (i_value < i_min ? i_min : (i_value > i_max ? i_max : i_value));
}

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
){
	const float inv = 1.f / (float) p_rows;
	for (unsigned j = 0; j < p_cols; j++) {
		float sum = 0.f;
		for (unsigned i = 0; i < p_rows; i++) {
			sum += io_group3d[j * p_rows + i];
		}

		o_baricenter[j] = sum * inv;

		for (unsigned i = 0; i < p_rows; i++) {
			io_group3d[j * p_rows + i] -= o_baricenter[j];
		}
	}
}

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
){
	float sigma = 0.f;

	for (unsigned c = 0; c < p_nChannels; c++) {
		//! Initialization
		float mean = 0.f;
		float std = 0.f;

		//! Compute the sum and the square sum
		for (unsigned n = 0; n < p_nSimP; n++) {
			for (unsigned k = 0; k < p_sP; k++) {
				const float value = i_Set[k + c * p_sP + n * p_sP * p_nChannels];
				mean += value;
				std  += value * value;
			}
		}

		//! Sample standard deviation (Bessel's correction)
		sigma += (std - mean * mean / (float) (p_sP * p_nSimP)) / (float) (p_sP * p_nSimP - 1);
	}

	return sigma / (float) p_nChannels;
}

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
){
    if (i_n == 1)
	 {
        o_a = 1;
        o_b = 1;
        return;
    }

    o_b = 2;
    while (i_n % o_b > 0) o_b++;
    o_a = i_n / o_b;

    if (o_b > o_a)
	 {
        o_a = o_b;
        o_b = i_n / o_a;
    }
}

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
){
    //! Open the file
    ofstream file;
    if (p_truncateFile) {
        file.open(p_pathName, ios::out | ios::trunc);
    }
    else {
        file.open(p_pathName, ios::out | ios::app);
    }

    //! Check if the file is open
    if (!file) {
        return EXIT_FAILURE;
    }

    //! Write measures in the file
    if (p_truncateFile) {
        file << "************" << endl;
        file << "-sigma = " << p_sigma << endl;
        file << "************" << endl;
    }
    file << "-PSNR" << p_app << " = " << p_psnr << endl;
    file << "-RMSE" << p_app << " = " << p_rmse << endl;
    file << "-GRPS" << p_app << " = " << p_grps << endl;
    file << "-TIME" << p_app << " = " << p_time << endl;
    file << "-CONS" << p_app << " = " << p_cons << endl;
    cout << endl;

    //! Close the file
    file.close();

    return EXIT_SUCCESS;
}

/**
 * @brief Check if a number is a power of 2
 **/
bool power_of_2(
    const unsigned n
){
    if (n == 0)
        return false;

    if (n == 1)
        return true;

    if (n % 2 == 0)
        return power_of_2(n / 2);
    else
        return false;
}

/**
 * @brief Look for the closest power of 2 number
 *
 * @param n: number
 *
 * @return the closest power of 2 lower or equal to n
 **/
int closest_power_of_2(
    const unsigned n
){
    unsigned r = 1;
    while (r * 2 <= n)
        r *= 2;

    return r;
}

/**
 * @brief Estimate sigma on each channel according to
 *        the choice of the color_space.
 *
 * @param sigma: estimated standard deviation of the noise;
 * @param sigma_Y : noise on the first channel;
 * @param sigma_U : (if chnls > 1) noise on the second channel;
 * @param sigma_V : (if chnls > 1) noise on the third channel;
 * @param chnls : number of channels of the image;
 * @param color_space : choice between OPP, YUV, YCbCr. If not
 *        then we assume that we're still in RGB space.
 *
 * @return EXIT_FAILURE if color_space has not expected
 *         type, otherwise return EXIT_SUCCESS.
 **/
int estimate_sigma(
    const float sigma
,   std::vector<float> &sigma_table
,   const unsigned chnls
,   const unsigned color_space
){
    if (chnls == 1)
        sigma_table[0] = sigma;
    else
    {
        if (color_space == YUV)
        {
            //! Y
            sigma_table[0] = sqrtf(0.299f * 0.299f + 0.587f * 0.587f + 0.114f * 0.114f) * sigma;
            //! U
            sigma_table[1] = sqrtf(0.14713f * 0.14713f + 0.28886f * 0.28886f + 0.436f * 0.436f) * sigma;
            //! V
            sigma_table[2] = sqrtf(0.615f * 0.615f + 0.51498f * 0.51498f + 0.10001f * 0.10001f) * sigma;
        }
        else if (color_space == YCBCR)
        {
            //! Y
            sigma_table[0] = sqrtf(0.299f * 0.299f + 0.587f * 0.587f + 0.114f * 0.114f) * sigma;
            //! U
            sigma_table[1] = sqrtf(0.169f * 0.169f + 0.331f * 0.331f + 0.500f * 0.500f) * sigma;
            //! V
            sigma_table[2] = sqrtf(0.500f * 0.500f + 0.419f * 0.419f + 0.081f * 0.081f) * sigma;
        }
        else if (color_space == OPP)
        {
            //! Y
            sigma_table[0] = sqrtf(0.333f * 0.333f + 0.333f * 0.333f + 0.333f * 0.333f) * sigma;
            //! U
            sigma_table[1] = sqrtf(0.5f * 0.5f + 0.0f * 0.0f + 0.5f * 0.5f) * sigma;
            //! V
            sigma_table[2] = sqrtf(0.25f * 0.25f + 0.5f * 0.5f + 0.25f * 0.25f) * sigma;
        }
        else if (color_space == RGB)
        {
            //! Y
            sigma_table[0] = sigma;
            //! U
            sigma_table[1] = sigma;
            //! V
            sigma_table[2] = sigma;
        }
        else
            return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


/**
 * @brief Initialize a 2D fftwf_plan with some parameters
 *
 * @param plan: fftwf_plan to allocate;
 * @param N: size of the patch to apply the 2D transform;
 * @param kind: forward or backward;
 * @param nb: number of 2D transform which will be processed.
 *
 * @return none.
 **/
void allocate_plan_2d(
    fftwf_plan* plan
,   const unsigned N
,   const fftwf_r2r_kind kind
,   const unsigned nb
){
    int            nb_table[2]   = {N, N};
    int            nembed[2]     = {N, N};
    fftwf_r2r_kind kind_table[2] = {kind, kind};

    float* vec = (float*) fftwf_malloc(N * N * nb * sizeof(float));
    (*plan) = fftwf_plan_many_r2r(2, nb_table, nb, vec, nembed, 1, N * N, vec,
                                  nembed, 1, N * N, kind_table, FFTW_ESTIMATE);

    fftwf_free(vec);
}

/**
 * @brief Initialize a 1D fftwf_plan with some parameters
 *
 * @param plan: fftwf_plan to allocate;
 * @param N: size of the vector to apply the 1D transform;
 * @param kind: forward or backward;
 * @param nb: number of 1D transform which will be processed.
 *
 * @return none.
 **/
void allocate_plan_1d(
    fftwf_plan* plan
,   const unsigned N
,   const fftwf_r2r_kind kind
,   const unsigned nb
){
    int nb_table[1] = {N};
    int nembed[1]   = {N * nb};
    fftwf_r2r_kind kind_table[1] = {kind};

    float* vec = (float*) fftwf_malloc(N * nb * sizeof(float));
    (*plan) = fftwf_plan_many_r2r(1, nb_table, nb, vec, nembed, 1, N, vec,
                                  nembed, 1, N, kind_table, FFTW_ESTIMATE);
    fftwf_free(vec);
}

/**
 * @brief Initialize a set of indices.
 *
 * @param ind_set: will contain the set of indices;
 * @param max_size: indices can't go over this size;
 * @param N : boundary;
 * @param step: step between two indices.
 *
 * @return none.
 **/
void ind_initialize(
    vector<unsigned> &ind_set
,   const unsigned beginning 
,   const unsigned end
,   const unsigned step
){
    ind_set.clear();
    unsigned ind = beginning;
    while (ind <= end)
    {
        ind_set.push_back(ind);
        ind += step;
    }
    if (ind_set.back() < end)
        ind_set.push_back(end);
}

void ind_initialize2(
    vector<unsigned> &ind_set
,   const unsigned max_size
,   const unsigned N
,   const unsigned step
){
    ind_set.clear();
    unsigned ind = N;
    while (ind < max_size - N)
    {
        ind_set.push_back(ind);
        ind += step;
    }
    if (ind_set.back() < max_size - N - 1)
        ind_set.push_back(max_size - N - 1);
}


/**
 * @brief For convenience. Estimate the size of the ind_set vector built
 *        with the function ind_initialize().
 *
 * @return size of ind_set vector built in ind_initialize().
 **/
unsigned ind_size(
    const unsigned beginning 
,   const unsigned end
,   const unsigned step
){
    unsigned ind = beginning;
    unsigned k = 0;
    while (ind <= end)
    {
        k++;
        ind += step;
    }
    if (ind - step < end)
        k++;

    return k;
}
