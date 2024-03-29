/*
 * Copyright (c) 2018, Thibaud Ehret <ehret.thibaud@gmail.com>
 * Copyright (c) 2018, Pablo Arias <pablo.arias@upf.edu>
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
 * @file vbm3d.cpp
 * @brief VBM3D denoising functions
 *
 * @author Thibaud Ehret <ehret.thibaud@gmail.com>
 **/

#define FLOAT_TRAJECTORIES

#include <iostream>
#include <algorithm>
#include <math.h>
#include <unordered_map>

#include "vbm3d.h"
#include "Utilities/Utilities.h"
#include "lib_transforms.h"

#define SQRT2     1.414213562373095
#define SQRT2_INV 0.7071067811865475
#define DCT       0
#define BIOR      1
#define HADAMARD  2
#define HAAR      4

#ifdef _OPENMP
#include <omp.h>
#endif

/* 
 * In order to reproduce the original VBM3D the DC coefficients are
 * not thresholded (DCTHRESH commented) but are filtered using Wiener
 * (DCWIENER uncommented), MTRICK activates undocumented tricks from
 * Marc Lebrun's implementation of BM3D available in IPOL
 * http://www.ipol.im/pub/art/2012/l-bm3d/
 */

//#define DCTHRESH
#define DCWIENER
//#define MTRICK

using namespace std;

bool ComparaisonFirst(pair<float,unsigned> pair1, pair<float,unsigned> pair2)
{
	return pair1.first < pair2.first;
}

/** ----------------- **/
/** - Main function - **/
/** ----------------- **/
/**
 * @brief run VBM3D process. Depending on whether OpenMP is used or not,
 *        and on the number of available threads, it divides the noisy
 *        video in sub_videos, to process them in parallel.
 *
 * @param sigma: value of assumed noise of the noisy video;
 * @param vid_noisy: noisy video;
 * @param vid_basic: will be the basic estimation after the 1st step
 * @param vid_denoised: will be the denoised final video;
 * @param useSD_h (resp. useSD_w): if true, use weight based
 *        on the standard variation of the 3D group for the
 *        first (resp. second) step, otherwise use the number
 *        of non-zero coefficients after Hard Thresholding
 *        (resp. the norm of Wiener coefficients);
 * @param tau_2D_hard (resp. tau_2D_wien): 2D transform to apply
 *        on every 3D group for the first (resp. second) part.
 *        Allowed values are DCT and BIOR;
 * @param color_space: Use a color transformation.
 *
 * @return EXIT_FAILURE if color_space has not expected
 *         type, otherwise return EXIT_SUCCESS.
 **/
int run_vbm3d(
    const float sigma
,   Video<float> &vid_noisy
,   Video<float> &fflow
,   Video<float> &bflow
,   Video<float> &vid_basic
,   Video<float> &vid_denoised
,   const Parameters& prms_1
,   const Parameters& prms_2
,   const bool color_space
){
	//! Check memory allocation
	if (vid_basic.sz != vid_noisy.sz)
		vid_basic.resize(vid_noisy.sz);
	if (vid_denoised.sz != vid_noisy.sz)
		vid_denoised.resize(vid_noisy.sz);

	const int chnls = vid_noisy.sz.channels;
 
	//! Transformation to OPP color space if necessary
	if(color_space)
	{
		printf("Transforming the color space\n");
		if(prms_1.k > 0)
			VideoUtils::transformColorSpace(vid_noisy, true);
		else
			VideoUtils::transformColorSpace(vid_basic, true);
	}

	//! Check if OpenMP is used or if number of cores of the computer is > 1
	unsigned nb_threads = 1;

#ifdef _OPENMP
	cout << "Open MP used" << endl;
	nb_threads = omp_get_max_threads();
	if(nb_threads > 32)
		nb_threads = 32;

	//! In case where the number of processors isn't a power of 2
	if (!power_of_2(nb_threads))
		nb_threads = closest_power_of_2(nb_threads);
#endif

	cout << endl << "Number of threads which will be used: " << nb_threads << endl;
#ifdef _OPENMP
	cout << " (real available cores: " << omp_get_max_threads() << ")" << endl;
#endif

	//! Allocate plan for FFTW library
	fftwf_plan plan_2d[nb_threads];
	fftwf_plan plan_2d_inv[nb_threads];
	fftwf_plan plan_1d[nb_threads];
	fftwf_plan plan_1d_inv[nb_threads];

	//! In the simple case
	if(nb_threads == 1)
	{
		Video<float> numerator, denominator;

		if(prms_1.k > 0)
		{
            const unsigned nb_cols = ind_size(0, vid_noisy.sz.width - prms_1.k, prms_1.p);
			//! Allocating Plan for FFTW process
			if (prms_1.T_2D == DCT)
			{
				allocate_plan_2d(&plan_2d[0], prms_1.k, FFTW_REDFT10,
						prms_1.N * chnls * prms_1.kt);
				allocate_plan_2d(&plan_2d_inv[0], prms_1.k, FFTW_REDFT01,
						prms_1.N * nb_cols * chnls * prms_1.kt);
			}
			allocate_plan_1d(&plan_1d[0], prms_1.kt, FFTW_REDFT10,
					prms_1.N * chnls * prms_1.k * prms_1.k);
			allocate_plan_1d(&plan_1d_inv[0], prms_1.kt, FFTW_REDFT01,
					prms_1.N * chnls * prms_1.k * prms_1.k);

			//! Denoising, 1st Step
			cout << "step 1..." << endl;
			vbm3d_1st_step(sigma, vid_noisy, fflow, bflow, prms_1,
					&plan_2d[0], &plan_2d_inv[0], &plan_1d[0], &plan_1d_inv[0],
					NULL, color_space, numerator, denominator);
			cout << "done." << endl;
			for (unsigned k = 0; k < vid_noisy.sz.whcf; k++)
				vid_basic(k) = (denominator(k) == 0) ? vid_noisy(k)
                                                : numerator(k) / denominator(k);
		}
		else
			cout << "skipping 1st step." << endl;


		if(prms_2.k > 0)
		{
			const unsigned nb_cols = ind_size(0, (vid_basic.sz.width - prms_2.k), prms_2.p);
			//! Allocating Plan for FFTW process
			if (prms_2.T_2D == DCT)
			{
				allocate_plan_2d(&plan_2d[0], prms_2.k, FFTW_REDFT10,
						prms_2.N * chnls * prms_2.kt);
				allocate_plan_2d(&plan_2d_inv[0], prms_2.k, FFTW_REDFT01,
						prms_2.N * nb_cols * chnls * prms_2.kt);
			}
			allocate_plan_1d(&plan_1d[0], prms_2.kt, FFTW_REDFT10,
					prms_2.N * chnls * prms_2.k * prms_2.k);
			allocate_plan_1d(&plan_1d_inv[0], prms_2.kt, FFTW_REDFT01,
					prms_2.N * chnls * prms_2.k * prms_2.k);

			//! Denoising, 2nd Step
			cout << "step 2..." << endl;
			vbm3d_2nd_step(sigma, vid_noisy, vid_basic, fflow, bflow,
					prms_2, &plan_2d[0], &plan_2d_inv[0], &plan_1d[0],
					&plan_1d_inv[0], NULL, color_space, numerator, denominator);
			cout << "done." << endl;
			for (unsigned k = 0; k < vid_noisy.sz.whcf; k++)
			{
				vid_denoised(k) = (denominator(k) == 0) ? vid_basic(k)
				                                        : numerator(k) / denominator(k);
			}
		}
		else
		{
			cout << "skipping 2nd step." << endl;
			for (unsigned k = 0; k < vid_noisy.sz.whcf; k++)
				vid_denoised(k) = vid_basic(k);
		}

	}
	//! If more than 1 threads are used
	else
	{
		//! Cut the video in nb_threads parts
		vector<Video<float> > numerator(nb_threads);
		vector<Video<float> > denominator(nb_threads);
		std::vector<VideoUtils::CropPosition > imCrops(nb_threads);
		VideoUtils::subDivideTight(vid_noisy, imCrops, 2 * prms_1.n, nb_threads);

		if(prms_1.k > 0)
		{
			//! Allocating Plan for FFTW process
			if (prms_1.T_2D == DCT)
				for (unsigned n = 0; n < nb_threads; n++)
				{
					const unsigned nb_cols = ind_size(
							(imCrops[n].origin_x == 0) ? 0 : prms_1.n,
							(imCrops[n].ending_x == imCrops[n].source_sz.width) ? 
								imCrops[n].ending_x - imCrops[n].origin_x - prms_1.k : 
								imCrops[n].ending_x - imCrops[n].origin_x - prms_1.k - prms_1.n,
							prms_1.p);
					allocate_plan_2d(&plan_2d[n], prms_1.k, FFTW_REDFT10,
							prms_1.N * chnls * prms_1.kt);
					allocate_plan_2d(&plan_2d_inv[n], prms_1.k, FFTW_REDFT01,
							prms_1.N * nb_cols * chnls * prms_1.kt);
				}
			for (unsigned n = 0; n < nb_threads; n++)
			{
				const unsigned nb_cols = ind_size(
						(imCrops[n].origin_x == 0) ? 0 : prms_1.n,
						(imCrops[n].ending_x == imCrops[n].source_sz.width) ?
							imCrops[n].ending_x - imCrops[n].origin_x - prms_1.k : 
							imCrops[n].ending_x - imCrops[n].origin_x - prms_1.k - prms_1.n,
							prms_1.p);
				allocate_plan_1d(&plan_1d[n], prms_1.kt, FFTW_REDFT10,
						prms_1.N * chnls * prms_1.k * prms_1.k);
				allocate_plan_1d(&plan_1d_inv[n], prms_1.kt, FFTW_REDFT01,
						prms_1.N * chnls * prms_1.k * prms_1.k);
			}

			//! denoising : 1st Step
			cout << "step 1..." << endl;;
#pragma omp parallel shared(\
		plan_2d, plan_2d_inv, plan_1d, plan_1d_inv, numerator, denominator, prms_1)
			{
#pragma omp for schedule(dynamic)
				for (unsigned n = 0; n < nb_threads; n++)
				{
					vbm3d_1st_step(sigma, vid_noisy, fflow, bflow, prms_1,
							&plan_2d[n], &plan_2d_inv[n], &plan_1d[n],
							&plan_1d_inv[n], &imCrops[n], color_space, numerator[n],
							denominator[n]);
				}
			}
			cout << "done." << endl;
			for (unsigned k = 0; k < vid_noisy.sz.whcf; k++)
			{
				float num = 0.;
				float den = 0.;
				for (unsigned n = 0; n < nb_threads; n++)
				{
					num += numerator[n](k);
					den += denominator[n](k);
				}
				if(den == 0)
					vid_basic(k) = vid_noisy(k);
				else
					vid_basic(k) = num / den;
			}
		}
		else
			cout << "skipping 1st step." << endl;

		if(prms_2.k > 0)
		{
			VideoUtils::subDivideTight(vid_basic, imCrops, 2 * prms_2.n, nb_threads);

			//! Allocating Plan for FFTW process
			if (prms_2.T_2D == DCT)
				for (unsigned n = 0; n < nb_threads; n++)
				{
					const unsigned nb_cols = ind_size(
							(imCrops[n].origin_x == 0) ? 0 : prms_2.n,
							(imCrops[n].ending_x == imCrops[n].source_sz.width) ?
								imCrops[n].ending_x - imCrops[n].origin_x - prms_2.k :
								imCrops[n].ending_x - imCrops[n].origin_x - prms_2.k - prms_2.n,
							prms_2.p);
					allocate_plan_2d(&plan_2d[n], prms_2.k, FFTW_REDFT10,
							prms_2.N * chnls * prms_2.kt);
					allocate_plan_2d(&plan_2d_inv[n], prms_2.k, FFTW_REDFT01,
							prms_2.N * nb_cols * chnls * prms_2.kt);
				}
			for (unsigned n = 0; n < nb_threads; n++)
			{
				const unsigned nb_cols = ind_size(
						(imCrops[n].origin_x == 0) ? 0 : prms_2.n,
						(imCrops[n].ending_x == imCrops[n].source_sz.width) ?
							imCrops[n].ending_x - imCrops[n].origin_x - prms_2.k :
							imCrops[n].ending_x - imCrops[n].origin_x - prms_2.k - prms_2.n,
						prms_2.p);
				allocate_plan_1d(&plan_1d[n], prms_2.kt, FFTW_REDFT10,
						prms_2.N * chnls * prms_2.k * prms_2.k);
				allocate_plan_1d(&plan_1d_inv[n], prms_2.kt, FFTW_REDFT01,
						prms_2.N * chnls * prms_2.k * prms_2.k);
			}

			//! Denoising: 2nd Step
			cout << "step 2..." << endl;
#pragma omp parallel shared(\
		plan_2d, plan_2d_inv, plan_1d, plan_1d_inv, numerator, denominator, prms_2)
			{
#pragma omp for schedule(dynamic)
				for (unsigned n = 0; n < nb_threads; n++)
				{
					vbm3d_2nd_step(sigma, vid_noisy, vid_basic, fflow, bflow,
							prms_2, &plan_2d[n], &plan_2d_inv[n], &plan_1d[n],
							&plan_1d_inv[n], &imCrops[n], color_space, numerator[n],
							denominator[n]);
				}
			}
			cout << "done." << endl;
			for (unsigned k = 0; k < vid_noisy.sz.whcf; k++)
			{
				float num = 0.;
				float den = 0.;
				for (unsigned n = 0; n < nb_threads; n++)
				{
					num += numerator[n](k);
					den += denominator[n](k);
				}
				if(den == 0)
					vid_denoised(k) = vid_noisy(k);
				else
					vid_denoised(k) = num / den;
			}
		}
		else
		{
			cout << "skipping 2nd step." << endl;
			for (unsigned k = 0; k < vid_noisy.sz.whcf; k++)
				vid_denoised(k) = vid_basic(k);
		}
	}

	//! Inverse color space transform to RGB if necessary
	if(color_space)
	{
		VideoUtils::transformColorSpace(vid_denoised,false);
		VideoUtils::transformColorSpace(vid_noisy, false);
		VideoUtils::transformColorSpace(vid_basic, false);
	}

	//! Free Memory
	if ((prms_1.k > 0 && prms_1.T_2D == DCT) || (prms_2.k > 0 && prms_2.T_2D == DCT))
		for(unsigned n = 0; n < nb_threads; n++)
		{
			fftwf_destroy_plan(plan_2d[n]);
			fftwf_destroy_plan(plan_2d_inv[n]);
			fftwf_destroy_plan(plan_1d[n]);
			fftwf_destroy_plan(plan_1d_inv[n]);
		}
	fftwf_cleanup();

	return EXIT_SUCCESS;
}

inline void patch_trajectory(
	int px
,	int py
,	int pt
,	const VideoSize &sz
,	const Video<float> &fflow
,	const Video<float> &bflow
,	const Parameters &prms
,	vector<vector<int> > &trajectory
){
	using std::min;
	using std::max;
	using std::round;
	const int k = prms.k;

#ifdef FLOAT_TRAJECTORIES
	// Initialize sub-pixel patch trajectories
	// (only needed if motion compensated patches are used)
	float fpx = px, fpy = py;
#endif

	for (int ht = 0; ht < prms.kt; ht++)
	{

		trajectory[ht][0] = px;
		trajectory[ht][1] = py;
		trajectory[ht][2] = pt + ht;

		if(prms.mc && ht < prms.kt - 1)
		{
#ifdef FLOAT_TRAJECTORIES
			// Next point in subpixel patch trajectory
			fpx += fflow(px + k/2, py + k/2, pt + ht, 0);
			fpy += fflow(px + k/2, py + k/2, pt + ht, 1);

			// Round to nearest integer
			px = (unsigned) min((int)sz.width  - k, max(0, (int)round(fpx)));
			py = (unsigned) min((int)sz.height - k, max(0, (int)round(fpy)));

#else // int trajectories
			int fflow_x = round(fflow(px + k/2, py + k/2, pt + ht, 0));
			int fflow_y = round(fflow(px + k/2, py + k/2, pt + ht, 1));
			px = (unsigned) min((int)sz.width  - k, max(0, (int)(px + fflow_x)));
			py = (unsigned) min((int)sz.height - k, max(0, (int)(py + fflow_y)));
#endif
		}
	}
}

/**
 * @brief Run the basic process of BM3D (1st step). The result
 *        is contained in vid_basic. The video has boundary, which
 *        are here only for block-matching and doesn't need to be
 *        denoised.
 *
 * @param sigma: value of assumed noise of the video to denoise;
 * @param vid_noisy: noisy video;
 * @param vid_basic: will contain the denoised video after the 1st step;
 * @param width, height, chnls : size of vid_noisy;
 * @param nHard: size of the boundary around vid_noisy;
 * @param useSD: if true, use weight based on the standard variation
 *        of the 3D group for the first step, otherwise use the number
 *        of non-zero coefficients after Hard-thresholding;
 * @param tau_2D: DCT or BIOR;
 * @param plan_2d_for_1, plan_2d_for_2, plan_2d_inv : for convenience. Used
 *        by fftw.
 *
 * @return none.
 **/
void vbm3d_1st_step(
    const float sigma
,   Video<float> const& vid_noisy
,   Video<float> &fflow
,   Video<float> &bflow
,   const Parameters& prms
,   fftwf_plan *  plan_2d
,   fftwf_plan *  plan_2d_inv
,   fftwf_plan *  plan_1d
,   fftwf_plan *  plan_1d_inv
,   VideoUtils::CropPosition* crop
,   const bool color_space
,   Video<float>& numerator
,   Video<float>& denominator
){
	VideoSize sz = vid_noisy.sz;

	//! Estimation of sigma on each channel
	const int chnls = sz.channels;
	vector<float> sigma_table(chnls);
	if (estimate_sigma(sigma, sigma_table, chnls, color_space) != EXIT_SUCCESS)
		return;

	//! Initialization for convenience
	vector<unsigned> row_ind(0);
	if(crop == NULL)
		ind_initialize(row_ind, 0, sz.height - prms.k, prms.p);
	else
		ind_initialize(row_ind, (crop->origin_y == 0) ? 0 : prms.n,
				(crop->ending_y == crop->source_sz.height) ?
					crop->ending_y - crop->origin_y - prms.k :
					crop->ending_y - crop->origin_y - prms.k - prms.n,
				prms.p);

	vector<unsigned> column_ind(0);
	if(crop == NULL)
		ind_initialize(column_ind, 0, sz.width - prms.k, prms.p);
	else
		ind_initialize(column_ind, (crop->origin_x == 0) ? 0 : prms.n,
				(crop->ending_x == crop->source_sz.width) ?
					crop->ending_x - crop->origin_x - prms.k :
					crop->ending_x - crop->origin_x - prms.k - prms.n,
				prms.p);

	const unsigned kHard_2 = prms.k * prms.k;
	const unsigned kHard_2_t = prms.k * prms.k * prms.kt;
	vector<float> group_3D_table(chnls * kHard_2_t * prms.N * column_ind.size());
	vector<float> wx_r_table;
	wx_r_table.reserve(chnls * column_ind.size());
	vector<float> hadamard_tmp(prms.N);

	//! Preprocessing (KaiserWindow, Threshold, DCT normalization, ...)
	vector<float> kaiser_window(kHard_2);
	vector<float> coef_norm(kHard_2);
	vector<float> coef_norm_inv(kHard_2);
	preProcess(kaiser_window, coef_norm, coef_norm_inv, prms.k);

	//! Preprocessing of Bior table
	vector<float> lpd, hpd, lpr, hpr;
	bior15_coef(lpd, hpd, lpr, hpr);

	//! Initialize aggregators
	denominator.resize(sz);
	numerator.resize(sz);
	for (unsigned k = 0; k < sz.whcf; k++)
	{
		numerator(k) = 0.f;
		denominator(k) = 0.f;
	}

	//! Precompute Block-Matching
	vector<vector<unsigned> > patch_table(column_ind.size(), std::vector<unsigned>(prms.N));
	vector<unsigned> size_patch_table(column_ind.size());
	vector<float> distances(prms.N);

	vector<float> table_2D(prms.N * chnls * kHard_2_t, 0.0f);

	vector<vector<int> > patch_traj(prms.kt, std::vector<int>(3));

	//! Loop on the frames
	for(unsigned t_r = 0; t_r < sz.frames - prms.kt + 1; ++t_r)
	{
		//! Loop on rows
		for (unsigned ind_i = 0; ind_i < row_ind.size(); ind_i++)
		{
			wx_r_table.clear();
			group_3D_table.clear();

			const unsigned i_r = row_ind[ind_i];

			//! Loop on columns
			for (unsigned ind_j = 0; ind_j < column_ind.size(); ind_j++)
			{
				//! Initialization
				const unsigned j_r = column_ind[ind_j];

				//! Number of similar patches
				unsigned pidx;
				if(crop == NULL) // video isn't a crop, nothing to change
					pidx = sz.index(j_r, i_r, t_r, 0);
				else // video is a crop: get position of the patch in original coordinates
					pidx = vid_noisy.getIndexSymmetric(crop->origin_x + j_r,
					                                   crop->origin_y + i_r,
					                                   crop->origin_t + t_r, 0);

				unsigned nSx_r = computeSimilarPatches(distances,
						patch_table[ind_j], pidx, vid_noisy, fflow, bflow, prms);
				size_patch_table[ind_j] = nSx_r;

				//! Update of table_2D
				if (prms.T_2D == DCT)
					dct_2d_process(table_2D, vid_noisy, patch_table[ind_j], plan_2d,
							coef_norm, fflow, prms);
				else if (prms.T_2D == BIOR)
					bior_2d_process(table_2D, vid_noisy, patch_table[ind_j], 
							lpd, hpd, fflow, prms);

				//! Build of the 3D group
				vector<float> group_3D(chnls * nSx_r * kHard_2_t, 0.0f);
				for(unsigned c = 0; c < chnls; c++)
					for (unsigned n = 0; n < nSx_r; n++)
						for (unsigned k = 0; k < kHard_2_t; k++)
							group_3D[n + k * nSx_r + c * kHard_2_t * nSx_r] =
								table_2D[k + n * kHard_2_t + c * kHard_2_t * prms.N];

                //! Transform along the temporal dimension
				if(prms.kt > 1 && prms.T_2D != BIOR)
					temporal_transform(group_3D, prms.k, prms.kt,
							chnls, nSx_r, prms.N, plan_1d);

				//! HT filtering of the 3D group
				vector<float> weight_table(chnls);
				if(prms.T_3D == HADAMARD)
					ht_filtering_hadamard(group_3D, hadamard_tmp, nSx_r, prms.k,
							prms.kt, chnls, sigma_table, prms.lambda3D, weight_table);
				else
					ht_filtering_haar(group_3D, hadamard_tmp, nSx_r, prms.k,
							prms.kt, chnls, sigma_table, prms.lambda3D, weight_table);

				//! Inverse transform along the temporal dimension
				if(prms.kt > 1 && prms.T_2D != BIOR)
					temporal_inv_transform(group_3D, prms.k, prms.kt, chnls, nSx_r,
							prms.N, plan_1d_inv);

				//! Save the 3D group. The DCT 2D inverse will be done after.
				for (unsigned c = 0; c < chnls; c++)
					for (unsigned n = 0; n < nSx_r; n++)
						for (unsigned k = 0; k < kHard_2_t; k++)
							group_3D_table.push_back(group_3D[n + k * nSx_r + c * kHard_2_t * nSx_r]);

				//! Save weighting
				for (unsigned c = 0; c < chnls; c++)
					wx_r_table.push_back(weight_table[c]);

			} //! loop on columns

			//!  Apply 2D inverse transform
			if (prms.T_2D == DCT)
				dct_2d_inv(group_3D_table, prms.k, prms.kt,
						prms.N * chnls * column_ind.size(), coef_norm_inv,
						plan_2d_inv);
			else if (prms.T_2D == BIOR)
				bior_2d_inv(group_3D_table, prms.k, prms.kt, lpr, hpr);


			//! Registration of the weighted estimation
			unsigned dec = 0;
			for (unsigned ind_j = 0; ind_j < column_ind.size(); ind_j++)
			{
				const unsigned nSx_r = size_patch_table[ind_j];
				for (unsigned c = 0; c < chnls; c++)
				{
					unsigned patch_j_r, patch_i_r, patch_t_r, patch_c_r;
					for (unsigned n = 0; n < nSx_r; n++)
					{
						sz.coords(patch_table[ind_j][n], patch_j_r,
						          patch_i_r, patch_t_r, patch_c_r);

						// compute patch trajectory
						patch_trajectory(patch_j_r, patch_i_r, patch_t_r, sz,
						                 fflow, bflow, prms, patch_traj);

						for (unsigned t = 0; t < prms.kt; t++)
						{
							// Coordinates of patch spatial slice at t
							patch_j_r = patch_traj[t][0];
							patch_i_r = patch_traj[t][1];

							for (unsigned p = 0; p < prms.k; p++)
							for (unsigned q = 0; q < prms.k; q++)
							{
								numerator(patch_j_r+q, patch_i_r+p, patch_t_r+t, c) +=
									kaiser_window[p * prms.k + q] * wx_r_table[c + ind_j * chnls]
									* group_3D_table[p * prms.k + q + t * kHard_2 + n * kHard_2_t +
									                 c * kHard_2_t * nSx_r + dec];
								denominator(patch_j_r+q, patch_i_r+p, patch_t_r+t, c) +=
									kaiser_window[p * prms.k + q] * wx_r_table[c + ind_j * chnls];
							}
						}
					}
				}

				dec += nSx_r * chnls * kHard_2_t;
			}

		} //! End of loop on i_r
	} //! End of loop on the frames
}

/**
 * @brief Run the final process of BM3D (2nd step). The result
 *        is contained in vid_denoised. The video has boundary, which
 *        are here only for block-matching and doesn't need to be
 *        denoised.
 *
 * @param sigma: value of assumed noise of the video to denoise;
 * @param vid_noisy: noisy video;
 * @param vid_basic: contains the denoised video after the 1st step;
 * @param vid_denoised: will contain the final estimate of the denoised
 *        video after the second step;
 * @param width, height, chnls : size of vid_noisy;
 * @param nWien: size of the boundary around vid_noisy;
 * @param useSD: if true, use weight based on the standard variation
 *        of the 3D group for the second step, otherwise use the norm
 *        of Wiener coefficients of the 3D group;
 * @param tau_2D: DCT or BIOR.
 *
 * @return none.
 **/
void vbm3d_2nd_step(
    const float sigma
,   Video<float> const& vid_noisy
,   Video<float> const& vid_basic
,   Video<float> &fflow
,   Video<float> &bflow
,   const Parameters& prms
,   fftwf_plan *  plan_2d
,   fftwf_plan *  plan_2d_inv
,   fftwf_plan *  plan_1d
,   fftwf_plan *  plan_1d_inv
,   VideoUtils::CropPosition* crop
,   const bool color_space
,   Video<float>& numerator
,   Video<float>& denominator
){
	VideoSize sz = vid_noisy.sz;

	//! Estimatation of sigma on each channel
	const int chnls = sz.channels;
	vector<float> sigma_table(chnls);
	if (estimate_sigma(sigma, sigma_table, chnls, color_space) != EXIT_SUCCESS)
		return;

	//! Initialization for convenience
	vector<unsigned> row_ind(0);
	if(crop == NULL)
		ind_initialize(row_ind, 0, (vid_basic.sz.height - prms.k), prms.p);
	else
		ind_initialize(row_ind, (crop->origin_y == 0) ? 0 : prms.n,
				(crop->ending_y == crop->source_sz.height) ?
					crop->ending_y - crop->origin_y - prms.k :
					crop->ending_y - crop->origin_y - prms.k - prms.n,
				prms.p);

	vector<unsigned> column_ind(0);
	if(crop == NULL)
		ind_initialize(column_ind, 0, (vid_basic.sz.width - prms.k), prms.p);
	else
		ind_initialize(column_ind, (crop->origin_x == 0) ? 0 : prms.n,
				(crop->ending_x == crop->source_sz.width) ?
					crop->ending_x - crop->origin_x - prms.k :
					crop->ending_x - crop->origin_x - prms.k - prms.n,
				prms.p);

	const unsigned kWien_2 = prms.k * prms.k;
	const unsigned kWien_2_t = kWien_2 * prms.kt;
	vector<float> group_3D_table(chnls * kWien_2_t * prms.N * column_ind.size());
	vector<float> wx_r_table;
	wx_r_table.reserve(chnls * column_ind.size());
	vector<float> tmp(prms.N);

	//! Preprocessing (KaiserWindow, Threshold, DCT normalization, ...)
	vector<float> kaiser_window(kWien_2);
	vector<float> coef_norm(kWien_2);
	vector<float> coef_norm_inv(kWien_2);
	preProcess(kaiser_window, coef_norm, coef_norm_inv, prms.k);

	//! Preprocessing of Bior table
	vector<float> lpd, hpd, lpr, hpr;
	bior15_coef(lpd, hpd, lpr, hpr);

	//! For aggregation part
	denominator.resize(sz);
	numerator.resize(sz);
	for (unsigned k = 0; k < sz.whcf; k++)
	{
		numerator(k) = 0.f;
		denominator(k) = 0.f;
	}

	//! Precompute Bloc-Matching
	vector<vector<unsigned> > patch_table(column_ind.size(), std::vector<unsigned>(prms.N));
	vector<unsigned> size_patch_table(column_ind.size());
	vector<float> distances(prms.N);

	//! DCT_table_2D[p * N + q + (i * width + j) * kWien_2 + c * (2 * ns + 1) * width * kWien_2]
	vector<float> table_2D_vid(prms.N * chnls * kWien_2_t, 0.0f);
	vector<float> table_2D_est(prms.N * chnls * kWien_2_t, 0.0f);

	vector<vector<int> > patch_traj(prms.kt, std::vector<int>(3));

	//§ Loop on the frames
	for(unsigned t_r = 0; t_r < sz.frames - prms.kt + 1; ++t_r)
	{
		//! Loop on i_r
		for (unsigned ind_i = 0; ind_i < row_ind.size(); ind_i++)
		{
			wx_r_table.clear();
			group_3D_table.clear();

			const unsigned i_r = row_ind[ind_i];

			//! Loop on j_r
			for (unsigned ind_j = 0; ind_j < column_ind.size(); ind_j++)
			{
				const unsigned j_r = column_ind[ind_j];

				//! Number of similar patches
				unsigned pidx;
				if(crop == NULL) // The video is not a cropped, nothing to change
					pidx = vid_basic.sz.index(j_r, i_r, t_r, 0);
				else // The video is a cropped, we need to find the original position of the patch
					pidx = vid_basic.getIndexSymmetric(crop->origin_x + j_r,
					                                   crop->origin_y + i_r,
					                                   crop->origin_t + t_r, 0);
				unsigned nSx_r = computeSimilarPatches(distances,
						patch_table[ind_j], pidx, vid_basic, fflow, bflow, prms);
				size_patch_table[ind_j] = nSx_r;

				//! Update of DCT_table_2D
				if (prms.T_2D == DCT)
				{
					dct_2d_process(table_2D_vid, vid_noisy, patch_table[ind_j], plan_2d,
							coef_norm, fflow, prms);
					dct_2d_process(table_2D_est, vid_basic, patch_table[ind_j], plan_2d,
							coef_norm, fflow, prms);
				}
				else if (prms.T_2D == BIOR)
				{
					bior_2d_process(table_2D_vid, vid_noisy, patch_table[ind_j],
							lpd, hpd, fflow, prms);
					bior_2d_process(table_2D_est, vid_basic, patch_table[ind_j],
							lpd, hpd, fflow, prms);
				}

				//! Build of the 3D group
				vector<float> group_3D_est(chnls * nSx_r * kWien_2_t, 0.0f);
				vector<float> group_3D_vid(chnls * nSx_r * kWien_2_t, 0.0f);
				for (unsigned c = 0; c < chnls; c++)
					for (unsigned n = 0; n < nSx_r; n++)
					{
						for (unsigned k = 0; k < kWien_2_t; k++)
						{
							group_3D_est[n + k * nSx_r + c * kWien_2_t * nSx_r] =
								table_2D_est[k + n * kWien_2_t + c * kWien_2_t * prms.N];
							group_3D_vid[n + k * nSx_r + c * kWien_2_t * nSx_r] =
								table_2D_vid[k + n * kWien_2_t + c * kWien_2_t * prms.N];
						}
					}

                //! Transform along the temporal dimension
				if(prms.kt > 1 && prms.T_2D != BIOR)
				{
					temporal_transform(group_3D_est, prms.k, prms.kt, chnls, nSx_r,
							prms.N, plan_1d);
					temporal_transform(group_3D_vid, prms.k, prms.kt, chnls, nSx_r,
							prms.N, plan_1d);
				}

				//! Wiener filtering of the 3D group
				vector<float> weight_table(chnls);
				if(prms.T_3D == HADAMARD)
					wiener_filtering_hadamard(group_3D_vid, group_3D_est, tmp,
							nSx_r, prms.k, prms.kt, chnls, sigma_table, weight_table);
				else
					wiener_filtering_haar(group_3D_vid, group_3D_est, tmp, nSx_r,
							prms.k, prms.kt, chnls, sigma_table, weight_table);

				//! Transform along the temporal dimension
				if(prms.kt > 1 && prms.T_2D != BIOR)
					temporal_inv_transform(group_3D_est, prms.k, prms.kt, chnls,
							nSx_r, prms.N, plan_1d_inv);

				//! Save the 3D group. The DCT 2D inverse will be done after.
				for (unsigned c = 0; c < chnls; c++)
				for (unsigned n = 0; n < nSx_r; n++)
				for (unsigned k = 0; k < kWien_2_t; k++)
					group_3D_table.push_back(group_3D_est[n + k * nSx_r + c * kWien_2_t * nSx_r]);

				//! Save weighting
				for (unsigned c = 0; c < chnls; c++)
					wx_r_table.push_back(weight_table[c]);

			} //! End of loop on j_r

			//!  Apply 2D dct inverse
			if (prms.T_2D == DCT)
				dct_2d_inv(group_3D_table, prms.k, prms.kt,
						prms.N * chnls * column_ind.size(),
						coef_norm_inv, plan_2d_inv);
			else if (prms.T_2D == BIOR)
				bior_2d_inv(group_3D_table, prms.k, prms.kt, lpr, hpr);

			//! Registration of the weighted estimation
			unsigned dec = 0;
			for (unsigned ind_j = 0; ind_j < column_ind.size(); ind_j++)
			{
				const unsigned nSx_r = size_patch_table[ind_j];
				for (unsigned c = 0; c < chnls; c++)
				{
					unsigned patch_j_r, patch_i_r, patch_t_r, patch_c_r;
					for (unsigned n = 0; n < nSx_r; n++)
					{
						sz.coords(patch_table[ind_j][n], patch_j_r,
						          patch_i_r, patch_t_r, patch_c_r);


//#ifdef FLOAT_TRAJECTORIES
//						// Initialize sub-pixel patch trajectory
//						// (only needed if motion compensated patches are used)
//						float fpatch_j_r = patch_j_r;
//						float fpatch_i_r = patch_i_r;
//#endif
						// compute patch trajectory
						patch_trajectory(patch_j_r, patch_i_r, patch_t_r, sz,
						                 fflow, bflow, prms, patch_traj);

						for (unsigned t = 0; t < prms.kt; t++)
						{
							// Coordinates of patch spatial slice at t
							patch_j_r = patch_traj[t][0];
							patch_i_r = patch_traj[t][1];

							for (unsigned p = 0; p < prms.k; p++)
							for (unsigned q = 0; q < prms.k; q++)
							{
								numerator(patch_j_r+q, patch_i_r+p, patch_t_r+t, c) +=
									kaiser_window[p * prms.k + q] * wx_r_table[c + ind_j * chnls] *
									group_3D_table[p * prms.k + q + t * kWien_2 + n * kWien_2_t +
									c * kWien_2_t * nSx_r + dec];
								denominator(patch_j_r+q, patch_i_r+p, patch_t_r+t, c) +=
									kaiser_window[p * prms.k + q] * wx_r_table[c + ind_j * chnls];
							}

//							if(prms.mc && t < prms.kt - 1)
//							{
//								using std::max;
//								using std::min;
//								using std::round;
//#ifdef FLOAT_TRAJECTORIES
//								// Next point in subpixel patch trajectory
//								fpatch_j_r += fflow(patch_j_r + prms.k/2, patch_i_r + prms.k/2, patch_t_r + t, 0);
//								fpatch_i_r += fflow(patch_j_r + prms.k/2, patch_i_r + prms.k/2, patch_t_r + t, 1);
//								// Round to nearest integer
//								patch_j_r = (unsigned)min(max((int)round(fpatch_j_r), 0), (int)(sz.width  - prms.k));
//								patch_i_r = (unsigned)min(max((int)round(fpatch_i_r), 0), (int)(sz.height - prms.k));
//#else // int trajectories
//								int fflow_j = round(fflow(patch_j_r + prms.k/2, patch_i_r + prms.k/2, patch_t_r + t, 0));
//								int fflow_i = round(fflow(patch_j_r + prms.k/2, patch_i_r + prms.k/2, patch_t_r + t, 1));
//								patch_j_r = (unsigned)min(max((int)(patch_j_r + fflow_j), 0), (int)(sz.width  - prms.k));
//								patch_i_r = (unsigned)min(max((int)(patch_i_r + fflow_i), 0), (int)(sz.height - prms.k));
//#endif
//							}
						}
					}
				}
				dec += nSx_r * chnls * kWien_2_t;
			}

		} //! End of loop on i_r
	}
}

inline float patchDistance(
	unsigned patch1
,	unsigned patch2
,	const Video<float>& vid
,	Video<float> &fflow
,	const Parameters &prms
,	vector<vector<int> > &traj_patch1
,	vector<vector<int> > &traj_patch2
){
	unsigned px, py, pt, pc;
	vid.sz.coords(patch1, px, py, pt, pc);

	unsigned qx, qy, qt, qc; 
	vid.sz.coords(patch2, qx, qy, qt, qc);

	using std::min;
	using std::max;
	using std::round;
	const int W = vid.sz.width;
	const int H = vid.sz.height;
	const int sPx = prms.k;
	const int sPt = prms.kt;

	float dist = 0.f, dif;
//#ifdef FLOAT_TRAJECTORIES
//	// Initialize sub-pixel patch trajectories
//	// (only needed if motion compensated patches are used)
//	float fpx = px, fpy = py;
//	float fqx = qx, fqy = qy;
//#endif

	// compute patch trajectories
	patch_trajectory(px, py, pt, vid.sz, fflow, Video<float>(), prms, traj_patch1);
	patch_trajectory(qx, qy, qt, vid.sz, fflow, Video<float>(), prms, traj_patch2);

	for (unsigned ht = 0; ht < sPt; ht++)
	{
		px = traj_patch1[ht][0], py = traj_patch1[ht][1];
		qx = traj_patch2[ht][0], qy = traj_patch2[ht][1];
		for (unsigned hc = 0; hc < vid.sz.channels; ++hc)
		for (unsigned hy = 0; hy < sPx; hy++)
		for (unsigned hx = 0; hx < sPx; hx++)
			dist += (dif = (vid(px + hx, py + hy, pt + ht, hc)
			              - vid(qx + hx, qy + hy, qt + ht, hc))) * dif;

//		if(mc && ht < sPt - 1)
//		{
//#ifdef FLOAT_TRAJECTORIES
//			// Next point in subpixel patch trajectory
//			fpx += fflow(px + sPx/2, py + sPx/2, pt + ht, 0);
//			fpy += fflow(px + sPx/2, py + sPx/2, pt + ht, 1);
//			// Round to nearest integer
//			px = (unsigned) min(W - sPx, max(0, (int)round(fpx)));
//			py = (unsigned) min(H - sPx, max(0, (int)round(fpy)));
//#else // int trajectories
//			int fflow_x = round(fflow(px + sPx/2, py + sPx/2, pt + ht, 0));
//			int fflow_y = round(fflow(px + sPx/2, py + sPx/2, pt + ht, 1));
//			px = (unsigned) min(W - sPx, max(0, (int)(px + fflow_x)));
//			py = (unsigned) min(H - sPx, max(0, (int)(py + fflow_y)));
//#endif
//		}
//
//#ifdef FLOAT_TRAJECTORIES
//		if(mc && ht < sPt - 1)
//		{
//			// Next point in subpixel patch trajectory
//			fqx += fflow(qx + sPx/2, qy + sPx/2, qt + ht, 0);
//			fqy += fflow(qx + sPx/2, qy + sPx/2, qt + ht, 1);
//			// Round to nearest integer
//			qx = (unsigned) min(W - sPx, max(0, (int)round(fqx)));
//			qy = (unsigned) min(H - sPx, max(0, (int)round(fqy)));
//#else // int trajectories
//			int fflow_x = round(fflow(qx + sPx/2, qy + sPx/2, qt + ht, 0))
//			int fflow_y = round(fflow(qx + sPx/2, qy + sPx/2, qt + ht, 1))
//			qx = (unsigned) min(W - sPx, max(0, (int)(qx + fflow_x)));
//			qy = (unsigned) min(H - sPx, max(0, (int)(qy + fflow_y)));
//#endif
//		}
	}

	return dist / (sPx * sPx * sPt * vid.sz.channels);
}

inline void localSearch(
	unsigned pidx
,	unsigned rpidx
,	unsigned s
,	const Video<float>& vid
,	std::unordered_map<unsigned, int>& alreadySeen
,	std::vector<std::pair<float, unsigned> >& bestPatches
,	Video<float> &fflow
,	const Parameters &prms
){
	int sWx = s;
	int sWy = s;
	int sPx = prms.k;
	int sPt = prms.kt;
	int Nb = prms.Nb;
	float d = prms.d;
	bool mc = prms.mc;

	//! Coordinates of the reference patch
	unsigned rpx, rpy, rpt, rpc;
	vid.sz.coords(rpidx, rpx, rpy, rpt, rpc);

	//! Coordinates of center of search box
	unsigned px, py, pt, pc;
	vid.sz.coords(pidx, px, py, pt, pc);

	unsigned rangex[2];
	unsigned rangey[2];

	vector<vector<int> > traj_patch1(sPt, vector<int>(3));
	vector<vector<int> > traj_patch2(sPt, vector<int>(3));

	rangex[0] = std::max(0, (int)px - (sWx-1)/2);
	rangey[0] = std::max(0, (int)py - (sWy-1)/2);

	rangex[1] = std::min((int)vid.sz.width  - sPx, (int)px + (sWx-1)/2);
	rangey[1] = std::min((int)vid.sz.height - sPx, (int)py + (sWy-1)/2);

	//! Redefine size of search range
	sWx = rangex[1] - rangex[0] + 1;
	sWy = rangey[1] - rangey[0] + 1;

	std::vector<std::pair<float, unsigned> > distance;
	distance.reserve(sWx*sWy);

	//! Compute distance between patches in search range
	for (unsigned qy = rangey[0], dy = 0; qy <= rangey[1]; qy++, dy++)
	for (unsigned qx = rangex[0], dx = 0; qx <= rangex[1]; qx++, dx++)
	{
		unsigned currentPatch = vid.sz.index(qx, qy, pt, 0);

		//! Save distance and corresponding patch index
		int seen = (alreadySeen[currentPatch]++);
		if(seen == 0)
		{
			float dist = patchDistance(rpidx, currentPatch, vid, fflow, prms,
			                           traj_patch1, traj_patch2)
			           - ((qx == rpx && qy == rpy) ? d : 0);
			distance.push_back(std::make_pair(dist, currentPatch));
		}
	}

	int nbCandidates = std::min(Nb, (int)distance.size());
	std::partial_sort(distance.begin(), distance.begin() + nbCandidates,
			distance.end(), comparaisonFirst);
	for(unsigned ix = 0; ix < nbCandidates; ++ix)
		bestPatches.push_back(distance[ix]);
}

int computeSimilarPatches(
	std::vector<float>& output
,	std::vector<unsigned>& index
,	unsigned pidx
,	const Video<float>& vid
,	Video<float> &fflow
,	Video<float> &bflow
,	const Parameters& prms
){
	std::vector<std::pair<float, unsigned> > bestPatches;
	bestPatches.reserve(prms.Nb*(2*prms.Nf+1));
	std::vector<unsigned> tempMatchesPre(prms.Nb);
	std::vector<unsigned> tempMatchesPost(prms.Nb);
	std::vector<std::pair<float, unsigned> > frameBestPatches;
	frameBestPatches.reserve(prms.Nb*prms.Nb);
	std::unordered_map<unsigned, int> alreadySeen;

	//! Coordinates of the reference patch
	unsigned rpx, rpy, rpt, rpc;
	vid.sz.coords(pidx, rpx, rpy, rpt, rpc);

	//! Coordinates of the entry points
	unsigned px, py, pt, pc;

	using std::min;
	using std::max;
	using std::round;
	int sPx = prms.k;
	int W = vid.sz.width;
	int H = vid.sz.height;

	//! Search in the current frame
	localSearch(pidx, pidx, prms.Ns, vid, alreadySeen, bestPatches, fflow, prms);
	for(unsigned ix = 0; ix < prms.Nb; ++ix)
	{
		tempMatchesPre[ix] = bestPatches[ix].second;
		tempMatchesPost[ix] = bestPatches[ix].second;
	}

	//! Search in the following frames (centered on the matches)
	int finalFrame = std::min(rpt + prms.Nf, vid.sz.frames - prms.kt) - rpt;	
	for(unsigned nextFrame = 0; nextFrame < finalFrame; ++nextFrame)
	{
		frameBestPatches.clear();
		for(unsigned currentTempMatch = 0; currentTempMatch < prms.Nb; ++currentTempMatch)
		{
			vid.sz.coords(tempMatchesPost[currentTempMatch], px, py, pt, pc);
			px = std::min(std::max((int)(px + std::round(fflow(px + prms.k/2,py + prms.k/2,pt,0))), 0), (int)(vid.sz.width - prms.k));
			py = std::min(std::max((int)(py + std::round(fflow(px + prms.k/2,py + prms.k/2,pt,1))), 0), (int)(vid.sz.height - prms.k));
			unsigned currentMatch = vid.sz.index(px, py, pt + 1, pc);
			localSearch(currentMatch, pidx, prms.Npr, vid, alreadySeen, frameBestPatches, fflow, prms);
		}

		int nbCandidates = std::min(prms.Nb, (unsigned)frameBestPatches.size());
		std::partial_sort(frameBestPatches.begin(), frameBestPatches.begin() + nbCandidates,
				frameBestPatches.end(), comparaisonFirst);
		for(unsigned ix = 0; ix < nbCandidates; ++ix)
		{
			tempMatchesPost[ix] = frameBestPatches[ix].second;
			bestPatches.push_back(frameBestPatches[ix]);
		}
	}

	//! Search in the previous frames (centered on the backward optical flow)
	finalFrame = rpt - std::max((int)(rpt - prms.Nf), 0);	
	for(unsigned nextFrame = 0; nextFrame < finalFrame; ++nextFrame)
	{
		frameBestPatches.clear();
		for(unsigned currentTempMatch = 0; currentTempMatch < prms.Nb; ++currentTempMatch)
		{
			vid.sz.coords(tempMatchesPre[currentTempMatch], px, py, pt, pc);
			px = std::min(std::max((int)(px + std::round(bflow(px + prms.k/2,py + prms.k/2,pt-1,0))), 0), (int)(vid.sz.width - prms.k));
			py = std::min(std::max((int)(py + std::round(bflow(px + prms.k/2,py + prms.k/2,pt-1,1))), 0), (int)(vid.sz.height - prms.k));
			unsigned currentMatch = vid.sz.index(px, py, pt - 1, pc);
			localSearch(currentMatch, pidx, prms.Npr, vid, alreadySeen, frameBestPatches, fflow, prms);
		}

		int nbCandidates = std::min(prms.Nb, (unsigned)frameBestPatches.size());
		std::partial_sort(frameBestPatches.begin(), frameBestPatches.begin() + nbCandidates,
				frameBestPatches.end(), comparaisonFirst);
		for(unsigned ix = 0; ix < nbCandidates; ++ix)
		{
			tempMatchesPre[ix] = frameBestPatches[ix].second;
			bestPatches.push_back(frameBestPatches[ix]);
		}
	}

	const unsigned nSimP = std::min(prms.N, (unsigned)bestPatches.size());

	std::partial_sort(bestPatches.begin(), bestPatches.begin() + nSimP,
	                  bestPatches.end(), comparaisonFirst);

	for (unsigned n = 0; n < nSimP; n++)
	{
		output[n] = bestPatches[n].first;
		index[n] = bestPatches[n].second;
	}

	unsigned ind_thresh = nSimP - 1;
	while((output[ind_thresh] > prms.tau) && (ind_thresh > 0))
		ind_thresh--;

	int candidates = closest_power_of_2(ind_thresh+1);

#ifdef MTRICK
	// Artificially adds a candidate when there's only the reference patch left
	if(candidates == 1)
	{
		candidates = 2;
		output[1] = output[0];
		index[1] = index[0];
	}
#endif

	return candidates;
}

/**
 * @brief Precompute a 2D DCT transform on all patches contained in
 *        a part of the video.
 *
 * @param DCT_table_2D : will contain the 2d DCT transform for all
 *        chosen patches;
 * @param vid : video on which the 2d DCT will be processed;
 * @param plan_1, plan_2 : for convenience. Used by fftw;
 * @param nHW : size of the boundary around vid;
 * @param width, height, chnls: size of vid;
 * @param kHW : size of patches (kHW x kHW);
 * @param i_r: current index of the reference patches;
 * @param step: space in pixels between two references patches;
 * @param coef_norm : normalization coefficients of the 2D DCT;
 * @param i_min (resp. i_max) : minimum (resp. maximum) value
 *        for i_r. In this case the whole 2d transform is applied
 *        on every patches. Otherwise the precomputed 2d DCT is re-used
 *        without processing it.
 **/
void dct_2d_process(
	vector<float> &DCT_table_2D
,	Video<float> const& vid
,	vector<unsigned> const& patch_table
,	fftwf_plan * plan
,	vector<float> const& coef_norm
,	Video<float> &fflow
,	const Parameters& prms
){
	//! Declarations
	const VideoSize sz = vid.sz;
	const unsigned kHW = prms.k;
	const unsigned ktHW = prms.kt;
	const unsigned kHW_2 = kHW * kHW;
	const unsigned kHW_2_t = kHW_2 * ktHW;
	const unsigned size = sz.channels * kHW_2_t * patch_table.size();

	vector<vector<int> > patch_traj(ktHW, vector<int>(3));

	using std::min;
	using std::max;
	using std::round;

	//! Allocating Memory
	float* vec = (float*) fftwf_malloc(size * sizeof(float));
	float* dct = (float*) fftwf_malloc(size * sizeof(float));

	for (unsigned c = 0; c < sz.channels; c++)
	{
		const unsigned dc_p = c * kHW_2_t * patch_table.size();
		for(unsigned n = 0; n < patch_table.size(); ++n)
		{
			unsigned i,j,t,tempc;
			sz.coords(patch_table[n], i, j, t, tempc);

//#ifdef FLOAT_TRAJECTORIES
//			// Initialize sub-pixel patch trajectories
//			// (only needed if motion compensated patches are used)
//			float fi = i, fj = j;
//#endif
			// compute patch trajectory
			patch_trajectory(i, j, t, sz, fflow, Video<float>(), prms, patch_traj);

			for (unsigned ht = 0; ht < ktHW; ht++)
			{
				i = patch_traj[ht][0];
				j = patch_traj[ht][1];
				for (unsigned p = 0; p < kHW; p++)
				for (unsigned q = 0; q < kHW; q++)
					vec[p * kHW + q + ht * kHW_2 + dc_p + n * kHW_2_t] = vid(i+q,j+p,t+ht,c);

//				if(mc && ht < ktHW - 1)
//				{
//#ifdef FLOAT_TRAJECTORIES
//					// Next point in subpixel patch trajectory
//					fi += fflow(i + kHW/2, j + kHW/2, t + ht, 0);
//					fj += fflow(i + kHW/2, j + kHW/2, t + ht, 1);
//					// Round to nearest integer
//					i = (unsigned)min(max((int)round(fi), 0), (int)(sz.width  - kHW));
//					j = (unsigned)min(max((int)round(fj), 0), (int)(sz.height - kHW));
//#else
//					int fflow_i = std::round(fflow(i + kHW/2, j + kHW/2, t + ht, 0));
//					int fflow_j = std::round(fflow(i + kHW/2, j + kHW/2, t + ht, 1));
//					i = (unsigned)min(max((int)(i + fflow_i), 0), (int)(sz.width  - kHW));
//					j = (unsigned)min(max((int)(j + fflow_j), 0), (int)(sz.height - kHW));
//#endif
//				}
			}
		}
	}

	//! Process of all DCTs
	fftwf_execute_r2r(*plan, vec, dct);
	fftwf_free(vec);

	//! Getting the result
	for (unsigned c = 0; c < vid.sz.channels; c++)
	{
		const unsigned dc_p = c * kHW_2_t * patch_table.size();
		for (unsigned n = 0; n < patch_table.size(); ++n)
		for (unsigned kt = 0; kt < ktHW; kt++)
		for (unsigned k = 0; k < kHW_2; k++)
			DCT_table_2D[dc_p + n * kHW_2_t + k + kt * kHW_2] =
				dct[dc_p + n * kHW_2_t + k + kt * kHW_2] * coef_norm[k];
	}
	fftwf_free(dct);
}

/**
 * @brief Precompute a 2D bior1.5 transform on all patches contained in
 *        a part of the video.
 *
 * @param bior_table_2D : will contain the 2d bior1.5 transform for all
 *        chosen patches;
 * @param vid : video on which the 2d transform will be processed;
 * @param nHW : size of the boundary around vid;
 * @param width, height, chnls: size of vid;
 * @param kHW : size of patches (kHW x kHW). MUST BE A POWER OF 2 !!!
 * @param i_r: current index of the reference patches;
 * @param step: space in pixels between two references patches;
 * @param i_min (resp. i_max) : minimum (resp. maximum) value
 *        for i_r. In this case the whole 2d transform is applied
 *        on every patches. Otherwise the precomputed 2d DCT is re-used
 *        without processing it;
 * @param lpd : low pass filter of the forward bior1.5 2d transform;
 * @param hpd : high pass filter of the forward bior1.5 2d transform.
 **/
void bior_2d_process(
	vector<float> &bior_table_2D
,	Video<float> const& vid
,	vector<unsigned> const& patch_table
,	vector<float> &lpd
,	vector<float> &hpd
,	Video<float> &fflow
,	const Parameters& prms
){
	//! Declarations
	const VideoSize sz = vid.sz;
	const unsigned kHW = prms.k;
	const unsigned ktHW = prms.kt;
	const unsigned kHW_2 = kHW * kHW;
	const unsigned kHW_2_t = kHW_2 * ktHW;

	vector<vector<int> > patch_traj(ktHW, vector<int>(3));

	//! If i_r == ns, then we have to process all Bior1.5 transforms
	for (unsigned c = 0; c < sz.channels; c++)
	{
		const unsigned dc_p = c * kHW_2_t * patch_table.size();
		for(unsigned n = 0; n < patch_table.size(); ++n)
		{
			unsigned i,j,t,tempc;
			sz.coords(patch_table[n], j, i ,t, tempc);

			patch_trajectory(j, i, t, sz, fflow, Video<float>(), prms, patch_traj);

//#ifdef FLOAT_TRAJECTORIES
//			// Initialize sub-pixel patch trajectories
//			// (only needed if motion compensated patches are used)
//			float fi = i, fj = j;
//#endif
			for(unsigned ht = 0; ht < ktHW; ++ht)
			{
				j = patch_traj[ht][0];
				i = patch_traj[ht][1];

				bior_2d_forward(vid, bior_table_2D, kHW, j, i, t+ht, c,
				                dc_p + n * kHW_2_t + ht * kHW_2, lpd, hpd);

//				if(prms.mc && ht < ktHW - 1)
//				{
//#ifdef FLOAT_TRAJECTORIES
//					// Next point in subpixel patch trajectory
//					fj += fflow(j + kHW/2, i + kHW/2, t + ht, 0);
//					fi += fflow(j + kHW/2, i + kHW/2, t + ht, 1);
//					// Round to nearest integer
//					j = (unsigned)min(max((int)round(fj), 0), (int)(sz.width  - kHW));
//					i = (unsigned)min(max((int)round(fi), 0), (int)(sz.height - kHW));
//#else
//					int fflow_j = std::round(fflow(j + kHW/2, i + kHW/2, t + ht, 0));
//					int fflow_i = std::round(fflow(j + kHW/2, i + kHW/2, t + ht, 1));
//					j = (unsigned)min(max((int)(j + fflow_j), 0), (int)(sz.width  - kHW));
//					i = (unsigned)min(max((int)(i + fflow_i), 0), (int)(sz.height - kHW));
//#endif
//				}
			}
		}
	}
}

/**
 * @brief HT filtering using Welsh-Hadamard transform (do only third
 *        dimension transform, Hard Thresholding and inverse transform).
 *
 * @param group_3D : contains the 3D block for a reference patch;
 * @param tmp: allocated vector used in Hadamard transform for convenience;
 * @param nSx_r : number of similar patches to a reference one;
 * @param kHW : size of patches (kHW x kHW);
 * @param chnls : number of channels of the video;
 * @param sigma_table : contains value of noise for each channel;
 * @param lambdaHard3D : value of thresholding;
 * @param weight_table: the weighting of this 3D group for each channel;
 * @param doWeight: if true process the weighting, do nothing
 *        otherwise.
 *
 * @return none.
 **/
void ht_filtering_hadamard(
    vector<float> &group_3D
,   vector<float> &tmp
,   const unsigned nSx_r
,   const unsigned kHard
,   const unsigned ktHard
,   const unsigned chnls
,   vector<float> const& sigma_table
,   const float lambdaHard3D
,   vector<float> &weight_table
){
	//! Declarations
	const unsigned kHard_2 = kHard * kHard;
	const unsigned kHard_2_t = kHard_2 * ktHard;
	for (unsigned c = 0; c < chnls; c++)
		weight_table[c] = 0.0f;
	const float coef_norm = sqrtf((float) nSx_r);
	const float coef = 1.0f / (float) nSx_r;

	//! Process the Welsh-Hadamard transform on the 3rd dimension
	for (unsigned n = 0; n < kHard_2_t * chnls; n++)
		hadamard_transform(group_3D, tmp, nSx_r, n * nSx_r);

	//! Hard Thresholding
	for (unsigned c = 0; c < chnls; c++)
	{
		const unsigned dc = c * nSx_r * kHard_2_t;
		const float T = lambdaHard3D * sigma_table[c] * coef_norm;
		for (unsigned k = 0; k < kHard_2_t * nSx_r; k++)
		{
#ifdef DCTHRESH
			if (fabs(group_3D[k + dc]) > T)
#else 
			if (k == 0 || fabs(group_3D[k + dc]) > T)
#endif
				weight_table[c]++;
			else
				group_3D[k + dc] = 0.0f;
		}
	}

	//! Process of the Welsh-Hadamard inverse transform
	for (unsigned n = 0; n < kHard_2_t * chnls; n++)
		hadamard_transform(group_3D, tmp, nSx_r, n * nSx_r);

	for (unsigned k = 0; k < group_3D.size(); k++)
		group_3D[k] *= coef;

	//! Weight for aggregation
	for (unsigned c = 0; c < chnls; c++)
		weight_table[c] = 1.0f / (float) (sigma_table[c] * sigma_table[c] * weight_table[c]);
}

/**
 * @brief HT filtering using Haar (do only third
 *        dimension transform, Hard Thresholding and inverse transform).
 *
 * @param group_3D : contains the 3D block for a reference patch;
 * @param tmp: allocated vector used in Hadamard transform for convenience;
 * @param nSx_r : number of similar patches to a reference one;
 * @param kHW : size of patches (kHW x kHW);
 * @param chnls : number of channels of the video;
 * @param sigma_table : contains value of noise for each channel;
 * @param lambdaHard3D : value of thresholding;
 * @param weight_table: the weighting of this 3D group for each channel;
 * @param doWeight: if true process the weighting, do nothing
 *        otherwise.
 *
 * @return none.
 **/
void ht_filtering_haar(
    vector<float> &group_3D
,   vector<float> &tmp
,   const unsigned nSx_r
,   const unsigned kHard
,   const unsigned ktHard
,   const unsigned chnls
,   vector<float> const& sigma_table
,   const float lambdaHard3D
,   vector<float> &weight_table
){
	//! Declarations
	const unsigned kHard_2 = kHard * kHard;
	const unsigned kHard_2_t = kHard_2 * ktHard;
	for (unsigned c = 0; c < chnls; c++)
		weight_table[c] = 0.0f;

	//! Process the Haar transform on the 3rd dimension
	for (unsigned n = 0; n < kHard_2_t * chnls; n++)
		haar_forward(group_3D, tmp, nSx_r, n * nSx_r);

	//! Hard Thresholding
	for (unsigned c = 0; c < chnls; c++)
	{
		const unsigned dc = c * nSx_r * kHard_2_t;
		const float T = lambdaHard3D * sigma_table[c];
		for (unsigned k = 0; k < kHard_2_t * nSx_r; k++)
		{
#ifdef DCTHRESH
			if (fabs(group_3D[k + dc]) > T)
#else
			if (k == 0 || fabs(group_3D[k + dc]) > T)
#endif
				weight_table[c]++;
			else
				group_3D[k + dc] = 0.0f;
		}
	}

	//! Process of the Haar inverse transform
	for (unsigned n = 0; n < kHard_2_t * chnls; n++)
		haar_inverse(group_3D, tmp, nSx_r, n * nSx_r);

	//! Weight for aggregation
	for (unsigned c = 0; c < chnls; c++)
		weight_table[c] = 1.f / (float) (sigma_table[c] * sigma_table[c] * weight_table[c]);
}

/**
 * @brief Wiener filtering using Hadamard transform.
 *
 * @param group_3D_vid : contains the 3D block built on vid_noisy;
 * @param group_3D_est : contains the 3D block built on vid_basic;
 * @param tmp: allocated vector used in hadamard transform for convenience;
 * @param nSx_r : number of similar patches to a reference one;
 * @param kWien : size of patches (kWien x kWien);
 * @param chnls : number of channels of the video;
 * @param sigma_table : contains value of noise for each channel;
 * @param weight_table: the weighting of this 3D group for each channel;
 * @param doWeight: if true process the weighting, do nothing
 *        otherwise.
 *
 * @return none.
 **/
void wiener_filtering_hadamard(
    vector<float> &group_3D_vid
,   vector<float> &group_3D_est
,   vector<float> &tmp
,   const unsigned nSx_r
,   const unsigned kWien
,   const unsigned ktWien
,   const unsigned chnls
,   vector<float> const& sigma_table
,   vector<float> &weight_table
){
	//! Declarations
	const unsigned kWien_2 = kWien * kWien;
	const unsigned kWien_2_t = kWien_2 * ktWien;
	const float coef = 1.0f / (float) nSx_r;

	for (unsigned c = 0; c < chnls; c++)
		weight_table[c] = 0.0f;

	//! Process the Welsh-Hadamard transform on the 3rd dimension
	for (unsigned n = 0; n < kWien_2_t * chnls; n++)
	{
		hadamard_transform(group_3D_vid, tmp, nSx_r, n * nSx_r);
		hadamard_transform(group_3D_est, tmp, nSx_r, n * nSx_r);
	}

	//! Wiener Filtering
	for (unsigned c = 0; c < chnls; c++)
	{
		const unsigned dc = c * nSx_r * kWien_2_t;
#ifdef DCWIENER
		for (unsigned k = 0; k < kWien_2_t * nSx_r; k++)
#else
		for (unsigned k = 1; k < kWien_2_t * nSx_r; k++)
#endif
		{
			float value = group_3D_est[dc + k] * group_3D_est[dc + k] * coef;
			value /= (value + sigma_table[c] * sigma_table[c]);
			group_3D_est[k + dc] = group_3D_vid[k + dc] * value * coef;
			weight_table[c] += (value*value);
		}
#ifndef DCWIENER
		// Add the weight corresponding to the DC components that was not thresholded
		weight_table[c] += nSx_r; 
#endif
	}

	//! Process of the Welsh-Hadamard inverse transform
	for (unsigned n = 0; n < kWien_2_t * chnls; n++)
		hadamard_transform(group_3D_est, tmp, nSx_r, n * nSx_r);

	//! Weight for aggregation
	for (unsigned c = 0; c < chnls; c++)
		weight_table[c] = (weight_table[c] > 0.0f ? 1.0f / (float)
				(sigma_table[c] * sigma_table[c] * weight_table[c]) : 1.0f);
}


/**
 * @brief Wiener filtering using Haar transform.
 *
 * @param group_3D_vid : contains the 3D block built on vid_noisy;
 * @param group_3D_est : contains the 3D block built on vid_basic;
 * @param tmp: allocated vector used in hadamard transform for convenience;
 * @param nSx_r : number of similar patches to a reference one;
 * @param kWien : size of patches (kWien x kWien);
 * @param chnls : number of channels of the video;
 * @param sigma_table : contains value of noise for each channel;
 * @param weight_table: the weighting of this 3D group for each channel;
 * @param doWeight: if true process the weighting, do nothing
 *        otherwise.
 *
 * @return none.
 **/
void wiener_filtering_haar(
    vector<float> &group_3D_vid
,   vector<float> &group_3D_est
,   vector<float> &tmp
,   const unsigned nSx_r
,   const unsigned kWien
,   const unsigned ktWien
,   const unsigned chnls
,   vector<float> const& sigma_table
,   vector<float> &weight_table
){
	//! Declarations
	const unsigned kWien_2 = kWien * kWien;
	const unsigned kWien_2_t = kWien_2 * ktWien;

	for (unsigned c = 0; c < chnls; c++)
		weight_table[c] = 0.0f;

	//! Process the Haar transform on the 3rd dimension
	for (unsigned n = 0; n < kWien_2_t * chnls; n++)
	{
		haar_forward(group_3D_vid, tmp, nSx_r, n * nSx_r);
		haar_forward(group_3D_est, tmp, nSx_r, n * nSx_r);
	}

	//! Wiener Filtering
	for (unsigned c = 0; c < chnls; c++)
	{
		const unsigned dc = c * nSx_r * kWien_2_t;
#ifdef DCWIENER
		for (unsigned k = 0; k < kWien_2_t * nSx_r; k++)
#else
		for (unsigned k = 1; k < kWien_2_t * nSx_r; k++)
#endif
		{
			float value = group_3D_est[dc + k] * group_3D_est[dc + k];
			value /= (value + sigma_table[c] * sigma_table[c]);
			group_3D_est[k + dc] = group_3D_vid[k + dc] * value;
			weight_table[c] += (value*value);
		}
#ifndef DCWIENER
		// Add the weight corresponding to the DC components that were not thresholded
		weight_table[c] += nSx_r; 
#endif
	}

	//! Process of the Welsh-Hadamard inverse transform
	for (unsigned n = 0; n < kWien_2_t * chnls; n++)
		haar_inverse(group_3D_est, tmp, nSx_r, n * nSx_r);

	//! Weight for aggregation
	for (unsigned c = 0; c < chnls; c++)
		weight_table[c] = 1.f / (float) (sigma_table[c] * sigma_table[c] * weight_table[c]);
}

/**
 * @brief Apply 2D dct inverse to a lot of patches.
 *
 * @param group_3D_table: contains a huge number of patches;
 * @param kHW : size of patch;
 * @param coef_norm_inv: contains normalization coefficients;
 * @param plan : for convenience. Used by fftw.
 *
 * @return none.
 **/
void dct_2d_inv(
    vector<float> &group_3D_table
,   const unsigned kHW
,   const unsigned ktHW
,   const unsigned N
,   vector<float> const& coef_norm_inv
,   fftwf_plan * plan
){
	//! Declarations
	const unsigned kHW_2 = kHW * kHW;
	const unsigned kHW_2_t = kHW_2 * ktHW;
	const unsigned size = kHW_2_t * N;
	const unsigned Ns   = group_3D_table.size() / kHW_2_t;

	//! Allocate Memory
	float* vec = (float*) fftwf_malloc(size * sizeof(float));
	float* dct = (float*) fftwf_malloc(size * sizeof(float));

	//! Normalization
	for (unsigned n = 0; n < Ns; n++)
		for (unsigned kt = 0; kt < ktHW; kt++)
			for (unsigned k = 0; k < kHW_2; k++)
				dct[k + n * kHW_2_t + kt * kHW_2] = group_3D_table[k + n * kHW_2_t + kt * kHW_2] * coef_norm_inv[k];

	//! 2D dct inverse
	fftwf_execute_r2r(*plan, dct, vec);
	fftwf_free(dct);

	//! Getting the result + normalization
	const float coef = 1.0f / (float)(kHW * 2);
	for (unsigned k = 0; k < group_3D_table.size(); k++)
		group_3D_table[k] = coef * vec[k];

	//! Free Memory
	fftwf_free(vec);
}

void bior_2d_inv(
    vector<float> &group_3D_table
,   const unsigned kHW
,   const unsigned ktHW
,   vector<float> const& lpr
,   vector<float> const& hpr
){
	//! Declarations
	const unsigned kHW_2 = kHW * kHW;
	const unsigned kHW_2_t = kHW_2 * ktHW;
	const unsigned N = group_3D_table.size() / kHW_2_t;

	//! Bior process
	for (unsigned n = 0; n < N; n++)
	for (unsigned kt = 0; kt < ktHW; kt++)
		bior_2d_inverse(group_3D_table, kHW, n * kHW_2_t + kt * kHW_2, lpr, hpr);
}

/** ----------------- **/
/** - Preprocessing - **/
/** ----------------- **/
/**
 * @brief Preprocess
 *
 * @param kaiser_window[kHW * kHW]: Will contain values of a Kaiser Window;
 * @param coef_norm: Will contain values used to normalize the 2D DCT;
 * @param coef_norm_inv: Will contain values used to normalize the 2D DCT;
 * @param bior1_5_for: will contain coefficients for the bior1.5 forward transform
 * @param bior1_5_inv: will contain coefficients for the bior1.5 inverse transform
 * @param kHW: size of patches (need to be 8 or 12).
 *
 * @return none.
 **/
void preProcess(
    vector<float> &kaiserWindow
,   vector<float> &coef_norm
,   vector<float> &coef_norm_inv
,   const unsigned kHW
){
	//! Kaiser Window coefficients
	if(kHW == 4)
	{
		//! First quarter of the matrix
		kaiserWindow[0 + kHW * 0] = 0.1924f; kaiserWindow[0 + kHW * 1] = 0.4055f;
		kaiserWindow[1 + kHW * 0] = 0.4055f; kaiserWindow[1 + kHW * 1] = 0.8544f;

		//! Completing the rest of the matrix by symmetry
		for(unsigned i = 0; i < kHW / 2; i++)
			for (unsigned j = kHW / 2; j < kHW; j++)
				kaiserWindow[i + kHW * j] = kaiserWindow[i + kHW * (kHW - j - 1)];

		for (unsigned i = kHW / 2; i < kHW; i++)
			for (unsigned j = 0; j < kHW; j++)
				kaiserWindow[i + kHW * j] = kaiserWindow[kHW - i - 1 + kHW * j];
	}
	else if (kHW == 6)
	{
		//! First quarter of the matrix
		kaiserWindow[0 + kHW * 0] = 0.1924f; kaiserWindow[0 + kHW * 1] = 0.3368f; kaiserWindow[0 + kHW * 2] = 0.4265f;
		kaiserWindow[1 + kHW * 0] = 0.3368f; kaiserWindow[1 + kHW * 1] = 0.5893f; kaiserWindow[1 + kHW * 2] = 0.7464f;
		kaiserWindow[2 + kHW * 0] = 0.4265f; kaiserWindow[2 + kHW * 1] = 0.7464f; kaiserWindow[2 + kHW * 2] = 0.9454f;

		//! Completing the rest of the matrix by symmetry
		for(unsigned i = 0; i < kHW / 2; i++)
			for (unsigned j = kHW / 2; j < kHW; j++)
				kaiserWindow[i + kHW * j] = kaiserWindow[i + kHW * (kHW - j - 1)];

		for (unsigned i = kHW / 2; i < kHW; i++)
			for (unsigned j = 0; j < kHW; j++)
				kaiserWindow[i + kHW * j] = kaiserWindow[kHW - i - 1 + kHW * j];
	}
	else if (kHW == 7)
	{
		//! First quarter of the matrix
		kaiserWindow[0 + kHW * 0] = 0.1924f; kaiserWindow[0 + kHW * 1] = 0.3151f; kaiserWindow[0 + kHW * 2] = 0.4055f; kaiserWindow[0 + kHW * 3] = 0.4387f;
		kaiserWindow[1 + kHW * 0] = 0.3151f; kaiserWindow[1 + kHW * 1] = 0.5161f; kaiserWindow[1 + kHW * 2] = 0.6640f; kaiserWindow[1 + kHW * 3] = 0.7184f;
		kaiserWindow[2 + kHW * 0] = 0.4055f; kaiserWindow[2 + kHW * 1] = 0.6640f; kaiserWindow[2 + kHW * 2] = 0.8544f; kaiserWindow[2 + kHW * 3] = 0.9243f;
		kaiserWindow[3 + kHW * 0] = 0.4387f; kaiserWindow[3 + kHW * 1] = 0.7184f; kaiserWindow[3 + kHW * 2] = 0.9243f; kaiserWindow[3 + kHW * 3] = 1.0000f; 

		//! Completing the rest of the matrix by symmetry
		for(unsigned i = 0; i <= kHW / 2; i++)
			for (unsigned j = kHW / 2 + 1; j < kHW; j++)
				kaiserWindow[i + kHW * j] = kaiserWindow[i + kHW * (kHW - j - 1)];

		for (unsigned i = kHW / 2 + 1; i < kHW; i++)
			for (unsigned j = 0; j < kHW; j++)
				kaiserWindow[i + kHW * j] = kaiserWindow[kHW - i - 1 + kHW * j];
	}
	else if (kHW == 8)
	{
		//! First quarter of the matrix
		kaiserWindow[0 + kHW * 0] = 0.1924f; kaiserWindow[0 + kHW * 1] = 0.2989f; kaiserWindow[0 + kHW * 2] = 0.3846f; kaiserWindow[0 + kHW * 3] = 0.4325f;
		kaiserWindow[1 + kHW * 0] = 0.2989f; kaiserWindow[1 + kHW * 1] = 0.4642f; kaiserWindow[1 + kHW * 2] = 0.5974f; kaiserWindow[1 + kHW * 3] = 0.6717f;
		kaiserWindow[2 + kHW * 0] = 0.3846f; kaiserWindow[2 + kHW * 1] = 0.5974f; kaiserWindow[2 + kHW * 2] = 0.7688f; kaiserWindow[2 + kHW * 3] = 0.8644f;
		kaiserWindow[3 + kHW * 0] = 0.4325f; kaiserWindow[3 + kHW * 1] = 0.6717f; kaiserWindow[3 + kHW * 2] = 0.8644f; kaiserWindow[3 + kHW * 3] = 0.9718f;

		//! Completing the rest of the matrix by symmetry
		for(unsigned i = 0; i < kHW / 2; i++)
			for (unsigned j = kHW / 2; j < kHW; j++)
				kaiserWindow[i + kHW * j] = kaiserWindow[i + kHW * (kHW - j - 1)];

		for (unsigned i = kHW / 2; i < kHW; i++)
			for (unsigned j = 0; j < kHW; j++)
				kaiserWindow[i + kHW * j] = kaiserWindow[kHW - i - 1 + kHW * j];
	}
	else if (kHW == 12)
	{
		//! First quarter of the matrix
		kaiserWindow[0 + kHW * 0] = 0.1924f; kaiserWindow[0 + kHW * 1] = 0.2615f; kaiserWindow[0 + kHW * 2] = 0.3251f; kaiserWindow[0 + kHW * 3] = 0.3782f;  kaiserWindow[0 + kHW * 4] = 0.4163f;  kaiserWindow[0 + kHW * 5] = 0.4362f;
		kaiserWindow[1 + kHW * 0] = 0.2615f; kaiserWindow[1 + kHW * 1] = 0.3554f; kaiserWindow[1 + kHW * 2] = 0.4419f; kaiserWindow[1 + kHW * 3] = 0.5139f;  kaiserWindow[1 + kHW * 4] = 0.5657f;  kaiserWindow[1 + kHW * 5] = 0.5927f;
		kaiserWindow[2 + kHW * 0] = 0.3251f; kaiserWindow[2 + kHW * 1] = 0.4419f; kaiserWindow[2 + kHW * 2] = 0.5494f; kaiserWindow[2 + kHW * 3] = 0.6390f;  kaiserWindow[2 + kHW * 4] = 0.7033f;  kaiserWindow[2 + kHW * 5] = 0.7369f;
		kaiserWindow[3 + kHW * 0] = 0.3782f; kaiserWindow[3 + kHW * 1] = 0.5139f; kaiserWindow[3 + kHW * 2] = 0.6390f; kaiserWindow[3 + kHW * 3] = 0.7433f;  kaiserWindow[3 + kHW * 4] = 0.8181f;  kaiserWindow[3 + kHW * 5] = 0.8572f;
		kaiserWindow[4 + kHW * 0] = 0.4163f; kaiserWindow[4 + kHW * 1] = 0.5657f; kaiserWindow[4 + kHW * 2] = 0.7033f; kaiserWindow[4 + kHW * 3] = 0.8181f;  kaiserWindow[4 + kHW * 4] = 0.9005f;  kaiserWindow[4 + kHW * 5] = 0.9435f;
		kaiserWindow[5 + kHW * 0] = 0.4362f; kaiserWindow[5 + kHW * 1] = 0.5927f; kaiserWindow[5 + kHW * 2] = 0.7369f; kaiserWindow[5 + kHW * 3] = 0.8572f;  kaiserWindow[5 + kHW * 4] = 0.9435f;  kaiserWindow[5 + kHW * 5] = 0.9885f;

		//! Completing the rest of the matrix by symmetry
		for(unsigned i = 0; i < kHW / 2; i++)
			for (unsigned j = kHW / 2; j < kHW; j++)
				kaiserWindow[i + kHW * j] = kaiserWindow[i + kHW * (kHW - j - 1)];

		for (unsigned i = kHW / 2; i < kHW; i++)
			for (unsigned j = 0; j < kHW; j++)
				kaiserWindow[i + kHW * j] = kaiserWindow[kHW - i - 1 + kHW * j];
	}
	else
		for (unsigned k = 0; k < kHW * kHW; k++)
			kaiserWindow[k] = 1.0f;

	//! Coefficient of normalization for DCT II and DCT II inverse
	const float coef = 0.5f / ((float) (kHW));
	for (unsigned i = 0; i < kHW; i++)
	for (unsigned j = 0; j < kHW; j++)
	{
		if (i == 0 && j == 0)
		{
			coef_norm    [i * kHW + j] = 0.5f * coef;
			coef_norm_inv[i * kHW + j] = 2.0f;
		}
		else if (i * j == 0)
		{
			coef_norm    [i * kHW + j] = SQRT2_INV * coef;
			coef_norm_inv[i * kHW + j] = SQRT2;
		}
		else
		{
			coef_norm    [i * kHW + j] = 1.0f * coef;
			coef_norm_inv[i * kHW + j] = 1.0f;
		}
	}
}

void temporal_transform(
    std::vector<float>& group_3D
,   const unsigned kHW
,   const unsigned ktHW
,   const unsigned chnls
,   const unsigned nSx_r
,   const unsigned N
,   fftwf_plan * plan
){
	//! Declarations
	const unsigned kHW_2 = kHW * kHW;
	const unsigned kHW_2_t = kHW_2 * ktHW;
    const float norm = 1./sqrt(2*ktHW);
    const unsigned size = chnls * kHW_2_t * N;

	//! Allocating Memory
	float* vec = (float*) fftwf_malloc(size * sizeof(float));
	float* dct = (float*) fftwf_malloc(size * sizeof(float));

    //! Load data
	for (unsigned c = 0; c < chnls; c++)
	{
		const unsigned dc_p = c * kHW_2_t * nSx_r;
		for (unsigned n = 0; n < nSx_r; ++n)
		for (unsigned k = 0; k < kHW_2; ++k)
		for (unsigned kt = 0; kt < ktHW; kt++)
			vec[k*ktHW + kt + dc_p + n * kHW_2_t] = group_3D[n + (k+kt*kHW_2)*nSx_r + c*kHW_2_t*nSx_r];
	}

	//! Process of all DCTs
	fftwf_execute_r2r(*plan, vec, dct);
	fftwf_free(vec);

	//! Getting the result
	for (unsigned c = 0; c < chnls; c++)
	{
		const unsigned dc_p = c * kHW_2_t * nSx_r;
		for (unsigned n = 0; n < nSx_r; ++n)
		for (unsigned k = 0; k < kHW_2; k++)
		for (unsigned kt = 0; kt < ktHW; kt++)
			group_3D[n + (k+kt*kHW_2)*nSx_r + c*kHW_2_t*nSx_r] =
				dct[k*ktHW + kt + dc_p + n * kHW_2_t] * norm * (kt==0?SQRT2_INV:1.);
	}
	fftwf_free(dct);
}

void temporal_inv_transform(
    std::vector<float>& group_3D
,   const unsigned kHW
,   const unsigned ktHW
,   const unsigned chnls
,   const unsigned nSx_r
,   const unsigned N
,   fftwf_plan * plan
){
	//! Declarations
	const unsigned kHW_2 = kHW * kHW;
	const unsigned kHW_2_t = kHW_2 * ktHW;
	const float norm = 1./sqrt(2*ktHW);
	const unsigned size = chnls * kHW_2_t * N;

	//! Allocating Memory
	float* vec = (float*) fftwf_malloc(size * sizeof(float));
	float* dct = (float*) fftwf_malloc(size * sizeof(float));

	//! Load data
	for (unsigned c = 0; c < chnls; c++)
	{
		const unsigned dc_p = c * kHW_2_t * nSx_r;
		for(unsigned n = 0; n < nSx_r; ++n)
		for(unsigned k = 0; k < kHW_2; ++k)
		for (unsigned kt = 0; kt < ktHW; kt++)
			dct[k*ktHW + kt + dc_p + n * kHW_2_t] = group_3D[n + (k+kt*kHW_2)*nSx_r + c*kHW_2_t*nSx_r] * (kt==0?SQRT2:1.);
	}

	//! Process of all DCTs
	fftwf_execute_r2r(*plan, dct, vec);
	fftwf_free(dct);

	//! Getting the result
	for (unsigned c = 0; c < chnls; c++)
	{
		const unsigned dc_p = c * kHW_2_t * nSx_r;
		for(unsigned n = 0; n < nSx_r; ++n)
		for (unsigned kt = 0; kt < ktHW; kt++)
		for (unsigned k = 0; k < kHW_2; k++)
			group_3D[n + (k+kt*kHW_2)*nSx_r + c*kHW_2_t*nSx_r] =
				vec[k*ktHW + kt + dc_p + n * kHW_2_t] * norm;
	}
	fftwf_free(vec);
}
