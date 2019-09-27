// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// Copyright (C) 2014, Nelson Monzón López <nmonzon@ctim.es>
// Copyright (C) 2014, Agustín Salgado de la Nuez <asalgado@dis.ulpgc.es>
// All rights reserved.

#ifndef GENERIC_TENSOR_H
#define GENERIC_TENSOR_H

#include <omp.h>

/**
 *
 * Compute the coefficients of the divergence term
 *
 */
void psi_divergence(
    float *psi1,      //coefficients of divergence
    float *psi2,      
    float *psi3,      
    float *psi4,      
    const float *psi, //robust functional
    const int nx,     //image width
    const int ny      //image height
)
{
    
    //calculate coefficients in the center body of the image
    #pragma omp parallel for
    for(int i = 1; i < ny-1; i++)
    {
	for(int j = 1; j < nx-1; j++)
	{
	    const int k = i * nx + j;

	    psi1[k] = 0.5   * (psi[k + 1] + psi[k]); 
	    psi2[k] = 0.5   * (psi[k - 1] + psi[k]);
	    psi3[k] = 0.5   * (psi[k +nx] + psi[k]);
	    psi4[k] = 0.5   * (psi[k -nx] + psi[k]);
	    
	}
    }      

    //calculate coefficients in the first and last rows
    #pragma omp parallel for
    for(int j = 1; j < nx-1; j++)
    {
	psi1[j] = 0.5   * (psi[j + 1] + psi[j]);
	psi2[j] = 0.5   * (psi[j - 1] + psi[j]);;
	psi3[j] = 0.5   * (psi[j +nx] + psi[j]);
	psi4[j] = 0; 

	const int k  = (ny-1)*nx + j;

	psi1[k] = 0.5   * (psi[k + 1] + psi[k]);;
	psi2[k] = 0.5   * (psi[k - 1] + psi[k]);
	psi3[k] = 0; 
	psi4[k] = 0.5   * (psi[k -nx] + psi[k]);
		    
    }    
    
    
      //calculate coefficients in the first and last columns
    #pragma omp parallel for
    for(int i = 1; i < ny-1; i++)
    {
	
	const int k  = i*nx;

	psi1[k] = 0.5   * (psi[k + 1] + psi[k]);
	psi2[k] = 0;  
	psi3[k] = 0.5   * (psi[k +nx] + psi[k]);
	psi4[k] = 0.5   * (psi[k -nx] + psi[k]);
	
	const int   j  = (i+1) * nx - 1;

	psi1[j] = 0; 
	psi2[j] = 0.5   * (psi[j - 1] + psi[j]);
	psi3[j] = 0.5   * (psi[j +nx] + psi[j]);
	psi4[j] = 0.5   * (psi[j -nx] + psi[j]);
	
    }
    
    
    //up-left corner (0,0) [0][0]
     psi1[0] = 0.5  * (psi[1] + psi[0]);
     psi3[0] = 0.5  * (psi[nx]  + psi[0]);
     psi2[0] = psi4[0] = 0;    

    //up-right corner (nx,0)  [0][nx-1]
     psi1[nx-1] = psi4[nx-1] = 0;
     psi2[nx-1] = 0.5   * (psi[nx-2] + psi[nx-1]);    
     psi3[nx-1] = 0.5   * (psi[2*nx-1] + psi[nx-1]);
      
    
    //bottom-left corner (0,ny) [ny-1][0]
     psi1[(ny-1)*nx] = 0.5 * (psi[(ny-1)*nx + 1] + psi[(ny-1)*nx]);
     psi4[(ny-1)*nx] = 0.5 * (psi[(ny-2)*nx] + psi[(ny-1) * nx]);
     psi2[(ny-1)*nx] = psi3[(ny-1)*nx] =  0;

    //bottom-right corner (nx,ny)   [ny-1][nx-1]
     psi2[ny*nx-1] = 0.5  * (psi[ny*nx - 2] + psi[ny*nx -1]);
     psi4[ny*nx-1] = 0.5  * (psi[ny*nx - 1 - nx] + psi[ny*nx -1]);
     psi1[ny*nx-1] = psi3[ny*nx-1] = 0;
}

/**
 *
 * Compute the divergence of the optical flow
 *
 */
void divergence(
    const float *u,    //component of optical flow
    const float *psi1, //coefficients of divergence
    const float *psi2, 
    const float *psi3, 
    const float *psi4, 
    const int nx,      //image width 
    const int ny,       //image height
    float *div      //computed divergence for u
)
{

    //calculate the divergence in the center body of the image
    #pragma omp parallel for
    for(int i = 1; i < ny-1; i++)
    {
	for(int j = 1; j < nx-1; j++)
	{
	    const int k = i * nx + j;

 	    div[k] = psi1[k] * (u[k + 1]    - u[k]) + psi2[k] * (u[k - 1]    - u[k]) + 
		       psi3[k] * (u[k + nx]   - u[k]) + psi4[k] * (u[k - nx]   - u[k]);
	}
    }

    //calculate the divergence in the first and last rows
    #pragma omp parallel for
    for(int j = 1; j < nx-1; j++)
    {

 	div[j] = psi1[j] * (u[j + 1]     - u[j]) + psi2[j] * (u[j - 1]     - u[j]) + 
 		   psi3[j] * (u[j + nx]    - u[j]);

	const int k = (ny-1)*nx + j;

 	div[k] = psi1[k] * (u[k + 1]     - u[k]) + psi2[k] * (u[k - 1]     - u[k]) + 
 		   psi4[k] * (u[k - nx]    - u[k]);
    }

    //calculate the divergence in the first and last columns
    #pragma omp parallel for
    for(int i = 1; i < ny-1; i++)
    {
	const int k = i*nx;

 	div[k] = psi1[k] * (u[k + 1]     - u[k]) + 
 		   psi3[k] * (u[k +nx]     - u[k]) + psi4[k] * (u[k - nx] - u[k]);
	
	const int j = (i+1) * nx - 1;

 	div[j] = psi2[j] * (u[j - 1]    - u[j]) + 
 		   psi3[j] * (u[j +nx]    - u[j]) + psi4[j] * (u[j -nx] - u[j]);
    }
    

    //up-left corner (0,0) [0][0]
    div[0] = psi1[0] * (u[1] - u[0]) + psi3[0] * (u[nx] - u[0]); 
    
    //up-right corner (nx,0) [0][nx]
    div[nx-1] = psi2[nx-1] * (u[nx - 2]   - u[nx-1]) + 
 		  psi3[nx-1] * (u[2*nx - 1] - u[nx-1]);
    
    //bottom-left corner (0,ny) [ny-1][0]
    div[(ny-1)*nx] = psi1[(ny-1)*nx] * (u[(ny-1) * nx + 1] - u[(ny-1) * nx]) + 
 		       psi4[(ny-1)*nx] * (u[(ny-2) * nx]     - u[(ny-1) * nx]);

    //bottom-right corner (nx,ny) [ny-1][nx-1]
     div[ny*nx-1] = psi2[ny*nx-1] * (u[ny*nx - 2]     - u[ny*nx -1]) + 
 		      psi4[ny*nx-1] * (u[ny*nx -1 -nx]  - u[ny*nx -1]);
}


#endif
