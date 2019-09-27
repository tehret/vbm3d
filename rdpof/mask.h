// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// Copyright (C) 2014, Nelson Monzón López  <nmonzon@ctim.es>
// All rights reserved.

#ifndef MASK_H
#define MASK_H

#include <omp.h>


/**
 *
 * Function to apply a 3x3 mask to an image
 *
 */
void
mask3x3 (const float *input,	//input image
	 float *output,		//output image
	 const int nx,		//image width
	 const int ny,		//image height
         const int nz,          // number of color channels in the image 
	 const float *mask	//mask to be applied
  )
{
  
  int nx_multichannel = nx * nz;
  
  for(int index_multichannel = 0; index_multichannel < nz; index_multichannel++){
  
	//apply the mask to the center body of the image
	
	#pragma omp parallel for
	for (int i = 1; i < ny - 1; i++){
		for (int j = 1; j < nx - 1; j++){
			
			int k = (i * nx + j) * nz + index_multichannel;
			double sum = 0;
			for (int l = 0; l < 3; l++){
				for (int m = 0; m < 3; m++){
					int p = ((i + l - 1) * nx + j + m - 1) * nz + index_multichannel;
					sum += input[p] * mask[l * 3 + m];
				}	
			}
			output[k] = sum;
		  }
    	}
      
   	//apply the mask to the first and last rows
        #pragma omp parallel for    	
	for (int j = 1; j < nx - 1; j++) {
	
    		int index = j * nz + index_multichannel;
		double sum = 0;
		
		sum += input[index - nz] * (mask[0] + mask[3]);
		sum += input[index] * (mask[1] + mask[4]);
		sum += input[index + nz] * (mask[2] + mask[5]);
        
		sum += input[nx_multichannel + j - nz] * mask[6];
		sum += input[nx_multichannel + j] * mask[7];
		sum += input[nx_multichannel + j + nz] * mask[8];

		output[j] = sum;
        
        	index = ((ny - 2) * nx + j) * nz + index_multichannel;
        
		sum = 0;
		sum += input[index - nz] * mask[0];
		sum += input[index] * mask[1];
		sum += input[index + nz] * mask[2];
        
       		index = ((ny - 1) * nx + j) * nz + index_multichannel;
        
		sum += input[index - nz] * (mask[6] + mask[3]);
		sum += input[index] * (mask[7] + mask[4]);
		sum += input[index + 1] * (mask[8] + mask[5]);

		output[index] = sum;
      	}
      
   	 //apply the mask to the first and last columns
	#pragma omp parallel for
        for (int i = 1; i < ny - 1; i++){
		
		int index = i * nx_multichannel + index_multichannel;
		
		double sum = 0;
		
		int index_row = (i - 1) * nx_multichannel + index_multichannel; 
		
		sum += input[index_row] * (mask[0] + mask[1]);
		sum += input[index_row + nz] * mask[2];

		sum += input[index] * (mask[3] + mask[4]);
		sum += input[index + nz] * mask[5];
        
        	index_row = (i + 1) * nx_multichannel + index_multichannel;
         
		sum += input[index_row] * (mask[6] + mask[7]);
		sum += input[index_row + nz] * mask[8];

		output[index] = sum;

		sum = 0;
		sum += input[index - 2 * nz] * mask[0];
		sum += input[index - nz] * (mask[1] + mask[2]);
       	 
        	index_row = (i + 1) * nx_multichannel + index_multichannel;
        
		sum += input[index_row - 2 * nz] * mask[3];
		sum += input[index_row - nz] * (mask[4] + mask[5]);
        
        	index_row = (i + 2) * nx_multichannel + index_multichannel;
        
		sum += input[index_row - 2 * nz] * mask[6];
		sum += input[index_row - nz] * (mask[7] + mask[8]);

		output[(i * nx + nx - 1) * nz + index_multichannel] = sum;

        }

    	//apply the mask to the four corners
	output[index_multichannel] = 
      		input[index_multichannel] * (mask[0] + mask[1] + mask[3] + mask[4]) + 
      		input[index_multichannel + nz] * (mask[2] + mask[5]) + 
      		input[nx_multichannel + index_multichannel] * (mask[6] + mask[7]) + 
      		input[nx_multichannel + index_multichannel + nz] * mask[8];

	output[nx_multichannel - nz + index_multichannel] = 
	  	input[(nx - 2) * nz + index_multichannel] * (mask[0] + mask[3]) + 
	  	input[(nx - 1) * nz + index_multichannel] * (mask[1] + mask[2] + mask[4] + mask[5]) + 
	  	input[(2 * nx - 2) * nz + index_multichannel] * mask[6] + 
	  	input[(2 * nx - 1) * nz + index_multichannel] * (mask[7] + mask[8]);

	output[(ny - 1) * nx_multichannel + index_multichannel] =
		input[(ny - 2) * nx_multichannel + index_multichannel] * (mask[0] + mask[1]) +
		input[((ny - 2) * nx + 1) * nz + index_multichannel] * mask[2] +
		input[(ny - 1) * nx_multichannel + index_multichannel] * (mask[3] + mask[4] + mask[6] + mask[7]) +
		input[((ny - 1) * nx + 1) * nz + index_multichannel] * (mask[5] + mask[8]);

	output[(ny * nx - 1) * nz + index_multichannel] =
		input[((ny - 1) * nx - 2) * nz + index_multichannel] * mask[0] +
		input[((ny - 1) * nx - 1) * nz + index_multichannel] * (mask[1] + mask[2]) +
		input[(ny * nx - 2) * nz + index_multichannel] * (mask[3] + mask[6]) +
		input[(ny * nx - 1) * nz + index_multichannel] * (mask[4] + mask[5] + mask[7] + mask[8]);
    	
     

  } // end loop for channels information
 
} // end mask3x3

/**
 *
 * Compute the second order X derivative
 *
 */
void
Dxx (const float *I,		//input image
     float *Ixx,		//oputput derivative
     const int nx,		//image width
     const int ny,		//image height
     const int nz               //number of color channels in the image  
  )
{
  //mask of second derivative
  float M[] = { 0., 0., 0.,
    1., -2., 1.,
    0., 0., 0.
  };

  //computing the second derivative
  mask3x3 (I, Ixx, nx, ny, nz, M);
}


/**
 *
 * Compute the second order Y derivative
 *
 */
void
Dyy (const float *I,		//input image
     float *Iyy,		//oputput derivative
     const int nx,		//image width
     const int ny,		//image height
     const int nz               //number of color channels in the image  
  )
{
  //mask of second derivative
  float M[] = { 0., 1., 0.,
    0., -2., 0.,
    0., 1., 0.
  };

  //computing the second derivative
  mask3x3 (I, Iyy, nx, ny, nz, M);
}


/**
 *
 * Compute the second order XY derivative
 *
 */
void
Dxy (const float *I,		//input image
     float *Ixy,		//oputput derivative
     const int nx,		//image width
     const int ny,		//image height
     const int nz               //number of color channels in the image  
  )
{
  //mask of second derivative
  float M[] = { 1. / 4., 0., -1. / 4.,
    0., 0., 0.,
    -1. / 4., 0., 1. / 4.
  };

  //computing the second derivative
  mask3x3 (I, Ixy, nx, ny, nz, M);
}


/**
 *
 * Compute the gradient with central differences
 *
 */
void
gradient (const float *input,	//input image
	  float *dx,		//computed x derivative
	  float *dy,		//computed y derivative
	  const int nx,		//image width
	  const int ny,		//image height
          const int nz          //number of color channels in the image 
  )
{
 
  const int nx_multichannel = nx * nz;

  for(int index_multichannel = 0; index_multichannel < nz; index_multichannel++){
	    
     //gradient in the center body of the image
     #pragma omp parallel for
     for (int i = 1; i < ny - 1; i++){
	for (int j = 1; j < nx - 1; j++){
		
		const int k = (i * nx + j) * nz + index_multichannel;
			
		dx[k] = 0.5 * (input[k + nz] - input[k - nz]);
		dy[k] = 0.5 * (input[k + nx_multichannel] - input[k - nx_multichannel]);
			
	  }
     }

    //gradient in the first and last rows
    #pragma omp parallel for    
    for (int j = 1; j < nx - 1; j++){
	    
	        const int index = j * nz + index_multichannel;
	    	
		dx[j] = 0.5 * (input[index + nz] - input[index - nz]);
		dy[j] = 0.5 * (input[index + nx_multichannel] - input[index]);

		const int k = ((ny - 1) * nx + j) * nz + index_multichannel;

		dx[k] = 0.5 * (input[k + nz] - input[k - nz]);
		dy[k] = 0.5 * (input[k] - input[k - nx_multichannel]);
     }

    //gradient in the first and last columns
    #pragma omp parallel for
    for (int i = 1; i < ny - 1; i++) {
		
		const int p = (i * nx_multichannel) + index_multichannel;
		
		dx[p] = 0.5 * (input[p + nz] - input[p]);
		dy[p] = 0.5 * (input[p + nx_multichannel] - input[p - nx_multichannel]);

		const int k = ((i + 1) * nx - 1) * nz + index_multichannel;

		dx[k] = 0.5 * (input[k] - input[k - nz]);
		dy[k] = 0.5 * (input[k + nx_multichannel] - input[k - nx_multichannel]);
    }
  
  
    //calculate the gradient in the corners
    dx[index_multichannel] = 0.5 * (input[index_multichannel + nz] - input[index_multichannel]);
    dy[index_multichannel] = 0.5 * (input[nx_multichannel + index_multichannel] - input[index_multichannel]);
    
    const int corner_up_right = (nx-1) * nz + index_multichannel;
    
    dx[corner_up_right] = 0.5 * (input[corner_up_right] - input[corner_up_right - nz]);
    dy[corner_up_right] = 0.5 * (input[(2 * nx_multichannel) + index_multichannel - nz] - input[corner_up_right]);
    
    const int corner_down_left = ((ny - 1) * nx) * nz + index_multichannel;
     
    dx[corner_down_left] = 0.5 * (input[corner_down_left + nz] - input[corner_down_left]);
    dy[corner_down_left] = 0.5 * (input[corner_down_left] - input[(ny - 2) * nx_multichannel + index_multichannel]) ;

    const int corner_down_right = ny * nx_multichannel - nz + index_multichannel;
    
    dx[corner_down_right] = 0.5 * (input[corner_down_right] - input[corner_down_right - nz]);
    dy[corner_down_right] = 0.5 * (input[corner_down_right] - input[(ny - 1) * nx_multichannel - nz + index_multichannel]);
  
  } // end loop for multi-channel
} // end gradient function


#endif
