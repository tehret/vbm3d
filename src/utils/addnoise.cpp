/*
 * Original work: Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * Modified work: Copyright (c) 2014, Pablo Arias <pariasm@gmail.com>
 * Modified work: Copyright (c) 2016, Thibaud Ehret <ehret.thibaud@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

#include <string>
#include <sstream>

#include "../Utilities/Utilities.h"
#include "../Utilities/LibVideoT.hpp"
#include "../Utilities/cmd_option.h"

using namespace std;

/**
 * @file   addnoise.cpp
 * @brief  Add noise to an image or video
 *
 * @author THIBAUD EHRET  <ehret.thibaud@gmail.com>
 **/

int main(int argc, char **argv)
{
	clo_usage("Add noise");
	clo_help(" NOTE: Input (<) and output (>) sequences are specified by their paths in printf format.\n");

	//! Paths to input/output sequences
	using std::string;
	const string  input_path = clo_option("-i", "", "< input sequence");
	const string  output_path = clo_option("-o", "noisy/noisy_%04d.tiff", "> output sequence");
	const int sigma = clo_option("-sigma", 20.f, "< Noise level (std. dev.)");
	const unsigned first_frame = clo_option("-f", 0, "< First frame of the video");
	const unsigned last_frame = clo_option("-l", 0, "< Last frame of the video");

	//! Declarations
	Video<float> original, noisy;

	//! Load input videos
	original.loadVideo(input_path, first_frame, last_frame, 1u);

	//! Add noise
	VideoUtils::addNoise(original, noisy, sigma, false);

	//! Save noisy video
	noisy.saveVideo(output_path, first_frame, 1u);
	return EXIT_SUCCESS;
}
