// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// Copyright (C) 2014, Agustín Salgado de la Nuez <asalgado@dis.ulpgc.es>
// Copyright (C) 2014, 2015 Nelson Monzón López <nmonzon@ctim.es>
// All rights reserved.

#ifndef SMOOTHNESS_H
#define SMOOTHNESS_H

#include <algorithm>
#include <vector>

#define EPSILON 0.001
#define XI 0.05
#define TAU 0.94 
#define BETA 0.001 

/**
  *
  * Compute the coefficients of the robust functional (smoothness term)
  *
**/
void psi_smooth(
    const float *ux, //gradient of x component of the optical flow
    const float *uy, //gradient of x component of the optical flow
    const float *vx, //gradient of y component of the optical flow
    const float *vy, //gradient of y component of the optical flow
    const float *expo, //exponential smoothing factor
    const int   size_flow,
    float       *psi       //output coefficients
)
{
    //compute 1/(sqrt(ux²+uy²+vx²+vy²+e²) in each pixel
    #pragma omp parallel for 
    for(int i = 0; i < size_flow; i++)
    {
	const float du  = expo[i] * ux[i] * ux[i] + expo[i] * uy[i] * uy[i];
	const float dv  = expo[i] * vx[i] * vx[i] + expo[i] * vy[i] * vy[i];
	const float normFlow  = du + dv;
	
	psi[i] = expo[i] / sqrt(normFlow + EPSILON * EPSILON);
			
    }
}

void max_gradients(
     const float *Ix, // Computed Image 1 derivative in x
     const float *Iy, // Computed Image 1 derivative in y
     const int   size, // Total image size (height * weight * nchannels)
     const int   nz, // nº channels
     float *maximum_gradients_per_pixel // vector with the maximum gradients per pixel
	
){
	int index_flow = 0;
	for(int index_image = 0; index_image < size; index_image+=nz){

		maximum_gradients_per_pixel[index_flow] = sqrt(Ix[index_image]*Ix[index_image] + Iy[index_image]*Iy[index_image]); 
				
		for(int index_multichannel = 1; index_multichannel < nz; index_multichannel++){
			  		
	  		const int real_index = index_image + index_multichannel;
			const float gradient_in_this_pixel = sqrt(Ix[real_index]*Ix[real_index] + Iy[real_index]*Iy[real_index]);		        
						
			if (maximum_gradients_per_pixel[index_flow] < gradient_in_this_pixel)
				maximum_gradients_per_pixel[index_flow] = gradient_in_this_pixel;	

        	}	
	 	index_flow++;
	}
        return;	
}

/**
  * Calculate the lambda optimum using the maximum gradient from all the multi-channel image	
  * It also return the maximum gradient per pixel 
**/
float lambda_optimum_using_maximum_gradient_per_pixel(
     const float *Ix, // Computed Image 1 derivative in x
     const float *Iy, // Computed Image 1 derivative in y
     const int   size, // Total image size (height * weight * nchannels)
     const int   size_flow, // Total flow size (height * weight)
     const int   nz, // nº channels
     const float alpha,  // smoothness weight
     float *lambda_per_pixel, // local lambda per pixel
     float *maximum_gradients_per_pixel // vector with the maximum gradients per pixel
)
{
	int index_flow = 0;
	
	for(int index_image = 0; index_image < size; index_image+=nz){

		maximum_gradients_per_pixel[index_flow] = sqrt(Ix[index_image]*Ix[index_image] + Iy[index_image]*Iy[index_image]); 
				
		for(int index_multichannel = 1; index_multichannel < nz; index_multichannel++){
			  		
	  		const int real_index = index_image + index_multichannel;
			const float gradient_in_this_pixel = sqrt(Ix[real_index]*Ix[real_index] + Iy[real_index]*Iy[real_index]);		        
						
			if (maximum_gradients_per_pixel[index_flow] < gradient_in_this_pixel)
				maximum_gradients_per_pixel[index_flow] = gradient_in_this_pixel;	

        	}	
		
	 	lambda_per_pixel[index_flow] = (-log(XI) + log(alpha)) / maximum_gradients_per_pixel[index_flow];
	 		
	 	index_flow++;
	}
	
	std::vector<float> gradients_ordered(size_flow);		
	#pragma omp parallel for 
        for (int index_flow = 0; index_flow < size_flow; index_flow++)
		gradients_ordered[index_flow] = maximum_gradients_per_pixel[index_flow];

	std::sort(gradients_ordered.begin(), gradients_ordered.end());

        const float c = -log(XI) + log(alpha);
	int pos_ref = TAU * size_flow;
	
	float lambda_pi = c / gradients_ordered[pos_ref-1];	
	
	while ((pos_ref < size_flow) && (c/2 > gradients_ordered[pos_ref-1]))
		pos_ref++;
	
	if (pos_ref == size_flow)
		lambda_pi = 0;
	else
		lambda_pi = (c / gradients_ordered[pos_ref-1]);

	return lambda_pi;
}

/**
 **  Calculate the exponential values.
 **
**/
void exponential_calculation(
		const float *Ix,         // Computed Image 1 derivative in x
		const float *Iy,         // Computed Image 1 derivative in y
		const int   size_flow,   // size of the flow field
		const int   size,        // size of the multi-channel image
		const int   nz,          // nº of image channels
		const float alpha,       // smoothness weight
		const float lambda,      // Coeffient for decreasing function
		const int   method_type, // (1 = DF, 2 = DF_BETA, 3 = DF_AUTO)
		float       *expo        // e^(lambda * DI)
		)
{

	float *maximum_gradients_per_pixel =  new float[size_flow];

	switch(method_type){

		case 1: case 2:{

			float beta = 0;
			if (method_type == 2) // For Exponential Beta Approximation
				beta = BETA;

			max_gradients(Ix, Iy, size, nz, maximum_gradients_per_pixel);

			#pragma omp parallel for
			for(int index_flow = 0; index_flow < size_flow; index_flow++){
				expo[index_flow] = exp(- lambda * maximum_gradients_per_pixel[index_flow]) + beta;
			}
			break;

		}
		case 3:{

			float *lambda_per_pixel = new float[size_flow];	
			float lambda_omega = lambda_optimum_using_maximum_gradient_per_pixel(Ix, Iy,
					size, size_flow, nz, alpha, lambda_per_pixel, maximum_gradients_per_pixel);

			#pragma omp parallel for
			for(int index_flow = 0; index_flow < size_flow; index_flow++){

				float lambda_pi = lambda_omega;
				if(lambda_omega > lambda_per_pixel[index_flow])
					lambda_pi = lambda_per_pixel[index_flow];

				expo[index_flow] = exp(- lambda_pi * maximum_gradients_per_pixel[index_flow]);

			}

			delete [] lambda_per_pixel;
			break;
		}
	} // End switch

	delete [] maximum_gradients_per_pixel;

}

#endif
