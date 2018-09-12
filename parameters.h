#ifndef PARAMETERS_H_INCLUDED
#define PARAMETERS_H_INCLUDED

/**
 * @brief Structures of parameters
 *
 **/

struct Parameters
{
	/// Type of the 2D tranform
	unsigned T_2D;
	/// Type of the 1D tranform 
	unsigned T_3D;
	/// Number of similar patches
	unsigned N;
	/// Number of frames forward (and backward) used during the search
	unsigned Nf;
	/// Size of the search region in the reference frame
	unsigned Ns;
	/// Size of the search region in the other frame
	unsigned Npr;
	/// Maximum number of matches kept for a frame
	unsigned Nb;
	/// Size of the patch (spatial)
	unsigned k;
	/// Size of the patch (temporal)
	unsigned kt;
	/// Step
	unsigned p;
	/// Correcting parameter in the distance computation
	float d;
	/// Threshold if it's a hard thresholding step
	float lambda3D;
	/// Distance threshold
	float tau;

	/// Border of the tile when using the multithreading
	int n = 16;
};

#endif
