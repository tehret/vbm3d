// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// Copyright (C) 2014, Agustín Salgado de la Nuez <asalgado@dis.ulpgc.es>
// Copyright (C) 2014, 2015 Nelson Monzón López <nmonzon@ctim.es>
// All rights reserved.

#ifndef ROBUST_EXPO_METHODS_H
#define ROBUST_EXPO_METHODS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <vector>
#include <iostream>
using namespace std;

#include "mask.h"
#include "zoom.h"
#include "bicubic_interpolation.h"
#include "generic_tensor.h"
#include "smoothness.h"

#define MAXITER 300
#define SOR_PARAMETER 1.9
#define GAUSSIAN_SIGMA 0.8


/**
  *
  * Compute the coefficients of the robust functional (data term)
  *
**/
void psi_data (
	const float *I1,	//first image
	const float *I2,	//second image
	const float *I2x,	//gradient of the second image
	const float *I2y,	//gradient of the second image
	const float *du,	//motion increment
	const float *dv,	//motion increment
	float *psip,		//output coefficients
	const int nx,		//image width
	const int ny,		//image height
	const int nz 		//image channels
)
{
	const int size = nx * ny;

	//compute 1/(sqrt(Sum(I2-I1+I2x*du+I2y*dv)²+e²) in each pixel
#pragma omp parallel for 
	for (int i = 0; i < size; i++)
	{
		float dI2 = 0;
		for (int k = 0; k < nz; k++){

			int   real_index = i*nz + k;
			const float dI  = I2[real_index] + I2x[real_index] * du[i] + I2y[real_index] * dv[i] - I1[real_index];

			dI2 += dI * dI;
		}

		psip[i] = (1. / sqrt (dI2 + EPSILON * EPSILON));

	}

}

/**
  *
  * Compute the coefficients of the robust functional (gradient term)
  *
**/

void psi_gradient (
	const float *I1x,	//gradient of the first image
	const float *I1y,		//gradient of the first image
	const float *I2x,		//gradient of the second image
	const float *I2y,		//gradient of the second image
	const float *I2xx,	//second derivatives of the second image
	const float *I2xy,	//second derivatives of the second image
	const float *I2yy,	//second derivatives of the second image
	const float *du,		//motion increment
	const float *dv,		//motion increment
	float *psip,		//output coefficients
	const int nx,		//image width
	const int ny,		//image height
	const int nz		// nº image channels
)
{
	const int size = nx * ny;

	//compute 1/(sqrt(|DI2-DI1=HI2*(du,dv)|²+e²) in each pixel
	#pragma omp parallel for
	for (int i = 0; i < size; i++)
	{
		float dI2 = 0;

		for (int k = 0; k < nz; k++){

			int   real_index = i*nz + k;

			const float dIx = I2x[real_index] + I2xx[real_index] * du[i] + I2xy[real_index] * dv[i] - I1x[real_index];
			const float dIy = I2y[real_index] + I2xy[real_index] * du[i] + I2yy[real_index] * dv[i] - I1y[real_index];

			dI2 += dIx * dIx + dIy * dIy;
		}

		psip[i] = (1. / sqrt (dI2 + EPSILON * EPSILON));

	}
}

/**
 * 
 *  SOR iteration in one position
 * 
 */
inline float sor_iteration(
    const float *Au,   //constant part of the numerator of u
    const float *Av,   //constant part of the numerator of vPriscila Estévez
    const float *Du,   //denominator of u
    const float *Dv,   //denominator of v
    const float *D,    //constant part of the numerator
    float       *du,   //x component of the motion increment
    float       *dv,   //y component of the motion increment
    const float alpha, //alpha smoothness parameter
    const float *psi1, //coefficients of the divergence
    const float *psi2, 
    const float *psi3,
    const float *psi4,
    const int   i,     //current row
    const int   i0,    //previous row
    const int   i1,    //following row
    const int   j,     //current column
    const int   j0,    //previous column
    const int   j1,    //following column
    const int   nx     //number of columns
)
{
	//set the SOR extrapolation parameter
	const float w = SOR_PARAMETER;

	//calculate the position in the array
	const int k = i * nx + j;

	//compute the divergence part of the numerator
	const float div_du = psi1[k] * du[k+j1] + psi2[k] * du[k-j0] +
	                     psi3[k] * du[k+i1] + psi4[k] * du[k-i0];

	const float div_dv = psi1[k] * dv[k+j1] + psi2[k] * dv[k-j0] +
	                     psi3[k] * dv[k+i1] + psi4[k] * dv[k-i0];

	const float duk = du[k];
	const float dvk = dv[k];

	//update the motion increment
	du[k] = (1.-w) * du[k] + w * (Au[k] - D[k] * dv[k] + alpha * div_du) / Du[k];
	dv[k] = (1.-w) * dv[k] + w * (Av[k] - D[k] * du[k] + alpha * div_dv) / Dv[k];

	//return the covergence error in this position
	return (du[k] - duk) * (du[k] - duk) + (dv[k] - dvk) * (dv[k] - dvk); 
}



/**
  *
  * Compute the optic flow with the matrix
  *
**/
void robust_expo_methods(
	const float *I1,          //first image
	const float *I2,          //second image
	float       *u,           //x component of the optical flow
	float       *v,           //y component of the optical flow
	const int    nx,          //image width
	const int    ny,          //image height
	const int    nz,          // number of color channels in the image
	const int    method_type, // choose the diffusion strategy
	const float  alpha,       // smoothness parameter
	const float  gamma,       // gradient term parameter
	const float  lambda,      // coefficient parameter for the decreasing function (if needed)
	const float  TOL,         // stopping criterion threshold
	const int    inner_iter,  // number of inner iterations
	const int    outer_iter,  // number of outer iterations
	const int    number_of_threads, // number of threads for the parallel code
	const bool   verbose      // switch on messages
)
{
	const int size_flow = nx * ny;
	const int size      = size_flow * nz;

	//allocate memory
	float *du    = new float[size_flow];
	float *dv    = new float[size_flow];

	float *ux    = new float[size_flow];
	float *uy    = new float[size_flow];
	float *vx    = new float[size_flow];
	float *vy    = new float[size_flow];

	float *I1x   = new float[size];
	float *I1y   = new float[size];
	float *I2x   = new float[size];
	float *I2y   = new float[size];
	float *I2w   = new float[size];
	float *I2wx  = new float[size];
	float *I2wy  = new float[size];
	float *I2xx  = new float[size];
	float *I2yy  = new float[size];
	float *I2xy  = new float[size];
	float *I2wxx = new float[size];
	float *I2wyy = new float[size];
	float *I2wxy = new float[size];

	float *div_u = new float[size_flow];
	float *div_v = new float[size_flow];
	float *div_d = new float[size_flow];

	float *Au    = new float[size_flow];
	float *Av    = new float[size_flow];
	float *Du    = new float[size_flow];
	float *Dv    = new float[size_flow];
	float *D     = new float[size_flow];
	float *expo  = new float[size_flow];

	float *psid  = new float[size_flow];
	float *psig  = new float[size_flow];
	float *psis  = new float[size_flow];
	float *psi1  = new float[size_flow];
	float *psi2  = new float[size_flow];
	float *psi3  = new float[size_flow];
	float *psi4  = new float[size_flow];

	//compute the gradient of the images
	gradient(I1, I1x, I1y, nx, ny, nz);
	gradient(I2, I2x, I2y, nx, ny, nz);

	//compute second order derivatives
	Dxx(I2, I2xx, nx, ny, nz);
	Dyy(I2, I2yy, nx, ny, nz);
	Dxy(I2, I2xy, nx, ny, nz);


	//compute the smoothness_coefficients including the robust function Phi for the smoothness term
	exponential_calculation (I1x, I1y, size_flow, size, nz, alpha, lambda, method_type, expo);

	//outer iterations loop
	for(int no = 0; no < outer_iter; no++){

		// Warp the second image and its derivatives
		bicubic_interpolation (I2, u, v, I2w, nx, ny, nz, true);
		bicubic_interpolation (I2x, u, v, I2wx, nx, ny, nz, true);
		bicubic_interpolation (I2y, u, v, I2wy, nx, ny, nz, true);
		bicubic_interpolation (I2xx, u, v, I2wxx, nx, ny, nz, true);
		bicubic_interpolation (I2xy, u, v, I2wxy, nx, ny, nz, true);
		bicubic_interpolation (I2yy, u, v, I2wyy, nx, ny, nz, true);

		// Compute the flow gradient
		gradient (u, ux, uy, nx, ny, 1);
		gradient (v, vx, vy, nx, ny, 1); 

		//compute robust function Phi for the smoothness term
		psi_smooth(ux, uy, vx, vy, expo, size_flow, psis);

		//compute coefficients of Phi functions in divergence
		psi_divergence (psi1, psi2, psi3, psi4, psis, nx, ny);

		//Calculate the divergence
		divergence (u, psi1, psi2, psi3, psi4, nx, ny, div_u);
		divergence (v, psi1, psi2, psi3, psi4, nx, ny, div_v);

#pragma omp parallel for
		for(int i = 0; i < size_flow; i++){

			//compute the coefficents of dw[i] in the smoothness term using gradient exponent
			div_d[i] = alpha * (psi1[i] + psi2[i] + psi3[i] + psi4[i]);

			//initialize the motion increment
			du[i] = dv[i] = 0;

		}

		//inner iterations loop
		for(int ni = 0; ni < inner_iter; ni++){

			//compute robust function Psi for the data and gradient terms
			psi_data (I1, I2w, I2wx, I2wy, du, dv, psid, nx, ny, nz);
			psi_gradient (I1x, I1y, I2wx, I2wy, I2wxx, I2wxy, I2wyy, du, dv, psig, nx, ny, nz);

			//store constant parts of the numerical scheme (equation (11))
			float BNu, dif, BNv, BDu, BDv, dx, dy, GNu, GNv, GDu, GDv, DI_Gradient, DI_Data;

			int index_flow = 0;

			for(int index_image = 0; index_image < size; index_image += nz){

				BNu = BNv = BDu = BDv = dx = dy = GNu = GNv = GDu = GDv = DI_Gradient = DI_Data = 0;

				for(int index_multichannel = 0; index_multichannel < nz; index_multichannel++){

					const int real_index = index_image + index_multichannel;

					//brightness constancy term
					dif  = I2w[real_index] - I1[real_index];
					BNu += dif * I2wx[real_index];
					BNv += dif * I2wy[real_index];
					BDu += I2wx[real_index] * I2wx[real_index];
					BDv += I2wy[real_index] * I2wy[real_index];
					DI_Data += (I2wy [real_index] * I2wx [real_index]);

					//gradient constancy term
					dx   = (I2wx[real_index] - I1x[real_index]);
					dy   = (I2wy[real_index] - I1y[real_index]);
					GNu += (dx * I2wxx[real_index] + dy * I2wxy[real_index]);
					GNv += (dx * I2wxy[real_index] + dy * I2wyy[real_index]);
					GDu += (I2wxx[real_index] * I2wxx[real_index] + I2wxy[real_index] * I2wxy[real_index]);
					GDv += (I2wyy[real_index] * I2wyy[real_index] + I2wxy[real_index] * I2wxy[real_index]);
					DI_Gradient  += (I2wxx[real_index] + I2wyy[real_index]) * I2wxy[real_index];

				}

				const float g = gamma * psig[index_flow];

				BNu  = -psid[index_flow] * BNu;
				BNv  = -psid[index_flow] * BNv;
				BDu  =  psid[index_flow] * BDu;
				BDv  =  psid[index_flow] * BDv;
				GNu  = -g * GNu;
				GNv  = -g * GNv;
				GDu  =  g * GDu;
				GDv  =  g * GDv;

				Au[index_flow] = BNu + GNu + alpha * div_u[index_flow];
				Av[index_flow] = BNv + GNv + alpha * div_v[index_flow];
				Du[index_flow] = BDu + GDu + div_d[index_flow];
				Dv[index_flow] = BDv + GDv + div_d[index_flow];
				D [index_flow] = psid[index_flow] * DI_Data + g * DI_Gradient;

				index_flow++;

			} // end image loop


			//sor iterations loop
			float error = 1000;
			int nsor = 0;

			while( error > TOL && nsor < MAXITER){

				error = 0;
				nsor++;

				//#pragma omp parallel for reduction(+:error)
#pragma omp parallel for reduction(+:error) num_threads((number_of_threads < (ny-3)) ? number_of_threads : (ny-3))
				//update the motion increment in the center of the images
				for(int i = 1; i < ny-1; i++){
					for(int j = 1; j < nx-1; j++)

						error += sor_iteration(
								Au, Av, Du, Dv, D, du, dv, alpha,
								psi1, psi2, psi3, psi4,
								i, nx, nx, j, 1, 1, nx
								);
				}

				//update the motion increment in the first and last rows
				for(int j = 1; j < nx-1; j++){

					error += sor_iteration(
							Au, Av, Du, Dv, D, du, dv, alpha,
							psi1, psi2, psi3, psi4,
							0, 0, nx, j, 1, 1, nx
							);

					error += sor_iteration(
							Au, Av, Du, Dv, D, du, dv, alpha,
							psi1, psi2, psi3, psi4,
							ny-1, nx, 0, j, 1, 1, nx
							);
				}

				//update the motion increment in the first and last columns
				for(int i = 1; i < ny-1; i++){

					error += sor_iteration(
							Au, Av, Du, Dv, D, du, dv, alpha,
							psi1, psi2, psi3, psi4,
							i, nx, nx, 0, 0, 1, nx
							);

					error += sor_iteration(
							Au, Av, Du, Dv, D, du, dv, alpha,
							psi1, psi2, psi3, psi4,
							i, nx, nx, nx-1, 1, 0, nx
							);
				}

				//process the top-left corner (0,0)
				error += sor_iteration(
						Au, Av, Du, Dv, D, du, dv, alpha,
						psi1, psi2, psi3, psi4,
						0, 0, nx, 0, 0, 1, nx
						);

				//process the top-right corner (0,nx-1)
				error += sor_iteration(
						Au, Av, Du, Dv, D, du, dv, alpha,
						psi1, psi2, psi3, psi4,
						0, 0, nx, nx-1, 1, 0, nx
						);

				//process the bottom-left corner (ny-1,0)
				error += sor_iteration(
						Au, Av, Du, Dv, D, du, dv, alpha,
						psi1, psi2, psi3, psi4,
						ny-1, nx, 0, 0, 0, 1, nx
						);

				//process the bottom-right corner (ny-1,nx-1)
				error += sor_iteration(
						Au, Av, Du, Dv, D, du, dv, alpha,
						psi1, psi2, psi3, psi4,
						ny-1, nx, 0, nx-1, 1, 0, nx
						);

				error = sqrt(error / size);

			}
			if(verbose)
				std::cout << "Iterations: " << nsor << " Error: " << error << std::endl;

		}

		//update the flow with the estimated motion increment
		for(int i = 0; i < size_flow; i++){
			u[i] += du[i];
			v[i] += dv[i];
		}

	}
	//delete allocated memory
	delete []du;
	delete []dv;

	delete []ux;
	delete []uy;
	delete []vx;
	delete []vy;

	delete []I1x;
	delete []I1y;
	delete []I2x;
	delete []I2y;
	delete []I2w;
	delete []I2wx;
	delete []I2wy;
	delete []I2xx;
	delete []I2yy;
	delete []I2xy;
	delete []I2wxx;
	delete []I2wyy;
	delete []I2wxy;

	delete []div_u;
	delete []div_v;
	delete []div_d;

	delete []Au;
	delete []Av;
	delete []Du;
	delete []Dv;
	delete []D;

	delete []expo;

	delete []psid;
	delete []psig;
	delete []psis;
	delete []psi1;
	delete []psi2;
	delete []psi3;
	delete []psi4;

}



/**
  *
  * Function to normalize the images between 0 and 255
  *
**/
void image_normalization (
	const float *I1,   //input image 1
	const float *I2,   //input image 2
	float *I1n,        //normalized output image 1
	float *I2n,        //normalized output image 2 
	int size,          //size of the image
	int nz             //number of color channels in the images
)
{

	float max, min;

	// Initial min and max values of the images
	for(int index_multichannel = 0; index_multichannel < nz; index_multichannel++){

		if (I1[index_multichannel] > I2[index_multichannel])
			max = I1[index_multichannel];
		else
			max = I2[index_multichannel];

		if (I1[index_multichannel] < I2[index_multichannel])
			min = I1[index_multichannel];
		else
			min = I2[index_multichannel];


		//compute the max and min values of the images
		for (int i = nz; i < size; i+=nz){
			int real_index = i + index_multichannel;
			if (I1[real_index] > max) max = I1[real_index];
			if (I1[real_index] < min) min = I1[real_index];
			if (I2[real_index] > max) max = I2[real_index];
			if (I2[real_index] < min) min = I2[real_index];
		}


		//compute the global max and min
		float den = max - min;

		if (den > 0){
			//normalize the images between 0 and 255
#pragma omp parallel for
			for (int index_image = 0; index_image < size; index_image +=nz)  {

				int real_index = index_image + index_multichannel;

				I1n[real_index] = 255.0 * (I1[real_index] - min) / den;
				I2n[real_index] = 255.0 * (I2[real_index] - min) / den;
			}
		}
		else{
			//copy the original data
#pragma omp parallel for
			for (int index_image = 0; index_image < size; index_image += nz)  {

				int real_index = index_image + index_multichannel;

				I1n[real_index] = I1[real_index];
				I2n[real_index] = I2[real_index];

			}
		}
	}
}



/**
  *
  *  Multiscale approach for computing the optical flow
  *
**/
void robust_expo_methods(
    const float *I1,          // first image
    const float *I2,          // second image
    float       *u,           // x component of the optical flow
    float       *v,           // y component of the optical flow
    const int    nxx,         // image width
    const int    nyy,         // image height
    const int    nzz,         // number of color channels in image
    const int    method_type, // choose the diffusion strategy
    const float  alpha,       // smoothness parameter
    const float  gamma,       // gradient term parameter
    const float  lambda,      // coefficient parameter for the decreasing function (if needed)
    const int    nscales,     // number of scales
    const int    fscale,      // finest scale
    const float  nu,          // downsampling factor
    const float  TOL,         // stopping criterion threshold
    const int    inner_iter,  // number of inner iterations
    const int    outer_iter,  // number of outer iterations
    const bool   verbose      // switch on messages
)
{
	int size = nxx * nyy * nzz;

	std::vector<float *> I1s(nscales);
	std::vector<float *> I2s(nscales);
	std::vector<float *> us (nscales);
	std::vector<float *> vs (nscales);

	std::vector<int> nx(nscales);
	std::vector<int> ny(nscales);

	I1s[0] = new float[size];
	I2s[0] = new float[size];

	//normalize the input images between 0 and 255
	image_normalization(I1, I2, I1s[0], I2s[0], size, nzz);

	//presmoothing the finest scale images
	gaussian(I1s[0], nxx, nyy, nzz, GAUSSIAN_SIGMA);
	gaussian(I2s[0], nxx, nyy, nzz, GAUSSIAN_SIGMA);

	us [0] = u;
	vs [0] = v;
	nx [0] = nxx;
	ny [0] = nyy;

	//create the scales
	for(int s = 1; s < nscales; s++){

		zoom_size(nx[s-1], ny[s-1], nx[s], ny[s], nu);
		const int sizes_flow = nx[s] * ny[s];
		const int sizes = sizes_flow * nzz;

		I1s[s] = new float[sizes];
		I2s[s] = new float[sizes];
		us[s]  = new float[sizes_flow];
		vs[s]  = new float[sizes_flow];

		//compute the zoom from the previous scale
		zoom_out(I1s[s-1], I1s[s], nx[s-1], ny[s-1], nzz, nu);
		zoom_out(I2s[s-1], I2s[s], nx[s-1], ny[s-1], nzz, nu);
	}

	//initialization of the optical flow at the coarsest scale
	for(int i = 0; i < nx[nscales-1] * ny[nscales-1]; i++)
		us[nscales-1][i] = vs[nscales-1][i] = 0.0;

	// readapt alpha for multichannel information
	int alpha_adapted_for_nchannels = alpha * nzz;

	int number_of_threads = 0;
#pragma omp parallel reduction(+:number_of_threads)
	number_of_threads += 1;

	//pyramidal approach for computing the optical flow
	for(int s = nscales-1; s >= 0; s--){

		if(verbose)
			std::cout << "Scale: " << s << std::endl;

		if (s >= fscale)
			robust_expo_methods(
					I1s[s], I2s[s], us[s], vs[s], nx[s], ny[s], nzz,
					method_type, alpha_adapted_for_nchannels, gamma, lambda,
					TOL, inner_iter, outer_iter, number_of_threads, verbose
					);

		//if it is not the finer scale, then upsample the optical flow and adapt it conveniently
		if(s){
			zoom_in_flow (us[s], us[s - 1], nx[s], ny[s], nx[s - 1], ny[s - 1]);
			zoom_in_flow (vs[s], vs[s - 1], nx[s], ny[s], nx[s - 1], ny[s - 1]);

			for(int i = 0; i < nx[s-1] * ny[s-1]; i++){
				us[s-1][i] *= 1.0 / nu;
				vs[s-1][i] *= 1.0 / nu;
			}
		}
	}

	//delete allocated memory
	delete []I1s[0];
	delete []I2s[0];

	for(int i = 1; i < nscales; i++){
		delete []I1s[i];
		delete []I2s[i];
		delete []us [i];
		delete []vs [i];
	}
}

#endif
