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

#ifndef PARAMETERS_H_INCLUDED
#define PARAMETERS_H_INCLUDED

/**
 * @brief Structures of parameters
 *
 **/

struct Parameters
{
	/// Size of the patch (spatial)
	unsigned k;
	/// Size of the patch (temporal)
	unsigned kt;
	/// Step size
	unsigned p;

	/// Maximum number of similar patches
	unsigned N;
	/// Maximum distance threshold
	float tau;
	/// Correcting parameter in the distance computation
	float d;

	/// Number of frames forward (and backward) used during the search
	unsigned Nf;
	/// Size of the search region in the reference frame
	unsigned Ns;
	/// Size of the search region in the other frame
	unsigned Npr;
	/// Maximum number of matches kept for a frame
	unsigned Nb;

	/// Type of the 2D tranform
	unsigned T_2D;
	/// Type of the 1D tranform 
	unsigned T_3D;
	/// Threshold if it's a hard thresholding step
	float lambda3D = 0.;

	/// Motion compensation for 3d patches
	bool mc = false;

	/// Border of the tile when using the multithreading
	int n = 16;
};

#endif
