#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string.h>

#include "vbm3d.h"
#include "Utilities/Utilities.h"
#include "Utilities/LibVideoT.hpp"
#include "Utilities/cmd_option.h"

#define YUV       0
#define YCBCR     1
#define OPP       2
#define RGB       3
#define DCT       4
#define BIOR      5
#define HADAMARD  6
#define HAAR      7
#define NONE      8

using namespace std;

void initializeParameters_1(
	Parameters& prms
,	const int k
,	const int kt
,	const int Nf
,	const int Ns
,	const int Npr
,	const int Nb
,	const int p
,	const int N
,	const int d
,	const float tau
,	const float lambda3D
,	const unsigned T_2D
,	const unsigned T_3D
,	const float sigma
){
	if(k < 0)
		prms.k = 8;
	else
		prms.k = k;

	if(kt < 0)
		prms.kt = 1;
	else
		prms.kt = kt;

	if(Nf < 0)
		prms.Nf = 4;
	else
		prms.Nf = Nf;
	if(Ns < 0)
		prms.Ns = 7;
	else
		prms.Ns = Ns;
	if(Npr < 0)
		prms.Npr = 5;
	else
		prms.Npr = Npr;
	if(Nb < 0)
		prms.Nb = 2;
	else
		prms.Nb = Nb;

	if(d < 0)
		prms.d = (7.*7.)/(255.*prms.k*prms.k);
	else
		prms.d = (d*d)/(255.*prms.k*prms.k);

	if(p < 0)
		prms.p = 6;
	else
		prms.p = p;

	if(N < 0)
		prms.N = 8;
	else
		prms.N = N;

	if(tau < 0)
		prms.tau = (sigma > 30) ? 4500/(255.*255.) : 3000/(255.*255.);
	else
		prms.tau = tau;

	if(lambda3D < 0)
		prms.lambda3D = 2.7f;
	else
		prms.lambda3D = lambda3D;

	if(T_2D == NONE)
		prms.T_2D = DCT;
	else
		prms.T_2D = T_2D;

	if(T_3D == NONE)
		prms.T_3D = HAAR;
	else
		prms.T_3D = T_3D;
}

void initializeParameters_2(
	Parameters& prms
,	const int k
,	const int kt
,	const int Nf
,	const int Ns
,	const int Npr
,	const int Nb
,	const int p
,	const int N
,	const int d
,	const float tau
,	const unsigned T_2D
,	const unsigned T_3D
,	const float sigma
){
	if(k < 0)
		prms.k = (sigma > 30) ? 8 : 7;
	else
		prms.k = k;

	if(kt < 0)
		prms.kt = 1;
	else
		prms.kt = kt;

	if(Nf < 0)
		prms.Nf = 4;
	else
		prms.Nf = Nf;

	if(Ns < 0)
		prms.Ns = 7;
	else
		prms.Ns = Ns;
	if(Npr < 0)
		prms.Npr = 5;
	else
		prms.Npr = Npr;
	if(Nb < 0)
		prms.Nb = 2;
	else
		prms.Nb = Nb;

	if(d < 0)
		prms.d = (3.*3.)/(255.*prms.k*prms.k);
	else
		prms.d = (d*d)/(255.*prms.k*prms.k);

	if(p < 0)
		prms.p = 4;
	else
		prms.p = p;

	if(N < 0)
		prms.N = 8;
	else
		prms.N = N;

	if(tau < 0)
		prms.tau = (sigma > 30) ? 3000/(255.*255.) : 1500/(255.*255.);
	else
		prms.tau = tau;

	if(T_2D == NONE)
		prms.T_2D = DCT;
	else
		prms.T_2D = T_2D;

	if(T_3D == NONE)
		prms.T_3D = HAAR;
	else
		prms.T_3D = T_3D;
}


/**
 * @file   main.cpp
 * @brief  Main executable file. Do not use lib_fftw to
 *         process DCT.
 *
 * @author MARC LEBRUN  <marc.lebrun@cmla.ens-cachan.fr>
 */


int main(int argc, char **argv)
{
    //! Check if there is the right call for the algorithm

	using std::string;
	const string  input_path = clo_option("-i"    , ""              , "< input sequence");
	const string  inbsc_path = clo_option("-b"    , ""              , "< input basic sequence");
	const string  noisy_path = clo_option("-nisy" , "nisy_%03d.tiff" , "> noisy sequence");
	const string  final_path = clo_option("-deno" , "deno_%03d.tiff" , "> denoised sequence");
	const string  basic_path = clo_option("-bsic" , "bsic_%03d.tiff" , "> basic denoised sequence");
	const string   diff_path = clo_option("-diff" , "diff_%03d.tiff" , "> difference sequence");
	const string   meas_path = clo_option("-meas" , "measure.txt"   , "> text file containing the measures");
#ifdef OPTICALFLOW
	const string   flow_path = clo_option("-flow" , "flow_%03d.flo" , "< optical flow ");
#endif

	const unsigned firstFrame = clo_option("-f", 0, "first frame");
	const unsigned lastFrame  = clo_option("-l", 0, "last frame");
	const unsigned frameStep  = clo_option("-s", 1, "frame step");

	//! General parameters
	const float fSigma = clo_option("-sigma", 0.f, "Add noise of standard deviation sigma");
	const bool addnoise  = (bool) clo_option("-add", true, "< add noise");
	const bool verbose  = (bool) clo_option("-verbose"     , true , "> verbose output");

	//! VBM3D parameters
	const int kHard = clo_option("-kHard", -1 , "< ");
	const int ktHard = clo_option("-ktHard", -1 , "< ");
	const int NfHard = clo_option("-NfHard", -1 , "< ");
	const int NsHard = clo_option("-NsHard", -1 , "< ");
	const int NprHard = clo_option("-NprHard", -1 , "< ");
	const int NbHard = clo_option("-NbHard", -1 , "< ");
	const int pHard = clo_option("-pHard", -1 , "< ");
	const int NHard = clo_option("-NHard", -1 , "< ");
	const int dHard = clo_option("-dHard", -1 , "< ");
	const float tauHard = clo_option("-tauHard", -1. , "< ");
	const int kWien = clo_option("-kWien", -1 , "< ");
	const int ktWien = clo_option("-ktWien", -1 , "< ");
	const int NfWien = clo_option("-NfWien", -1 , "< ");
	const int NsWien = clo_option("-NsWien", -1 , "< ");
	const int NprWien = clo_option("-NprWien", -1 , "< ");
	const int NbWien = clo_option("-NbWien", -1 , "< ");
	const int pWien = clo_option("-pWien", -1 , "< ");
	const int NWien = clo_option("-NWien", -1 , "< ");
	const int dWien = clo_option("-dWien", -1 , "< ");
	const float tauWien = clo_option("-tauWien", -1. , "< ");

	const float lambda3D = clo_option("-lambda3d", -1. , "< ");
	const unsigned color_space  =  (unsigned) clo_option("-color", 0 , "< color space");

	//! Check inputs
	if (input_path == "")
	{
		fprintf(stderr, "%s: no input images.\nTry `%s --help' for more information.\n",
				argv[0], argv[0]);
		return EXIT_FAILURE;
	}

	//! Variables initialization
	const unsigned T_2D_hard  = (unsigned) clo_option("-T2dh", NONE , "< tau_2D_hard");
	if (T_2D_hard != NONE && T_2D_hard != DCT && T_2D_hard != BIOR)
	{
		cout << "T_2d_hard is not known. Choice is :" << endl;
		cout << " -dct (" << DCT << ")" << endl;
		cout << " -bior (" << BIOR << ")" << endl;
		return EXIT_FAILURE;
	}

	const unsigned T_2D_wien  = (unsigned) clo_option("-T2dw", NONE , "< tau_2D_wien");
	if (T_2D_wien != NONE && T_2D_wien != DCT && T_2D_wien != BIOR)
	{
		cout << "T_2d_wien is not known. Choice is :" << endl;
		cout << " -dct (" << DCT << ")" << endl;
		cout << " -bior (" << BIOR << ")" << endl;
		return EXIT_FAILURE;
	};

	const unsigned T_3D_hard  = (unsigned) clo_option("-T3dh", NONE , "< tau_3D_hard");
	if (T_3D_hard != NONE && T_3D_hard != HAAR && T_3D_hard != HADAMARD)
	{
		cout << "T_3d_hard is not known. Choice is :" << endl;
		cout << " -haar (" << HAAR << ")" << endl;
		cout << " -hadamard (" << HADAMARD << ")" << endl;
		return EXIT_FAILURE;
	}

	const unsigned T_3D_wien  = (unsigned) clo_option("-T3dw", NONE , "< tau_3D_wien");
	if (T_3D_wien != NONE && T_3D_wien != HAAR && T_3D_wien != HADAMARD)
	{
		cout << "T_3d_wien is not known. Choice is :" << endl;
		cout << " -haar (" << HAAR << ")" << endl;
		cout << " -hadamard (" << HADAMARD << ")" << endl;
		return EXIT_FAILURE;
	};

	Parameters prms_1;
	Parameters prms_2;

	initializeParameters_1(prms_1, kHard, ktHard, NfHard, NsHard, NprHard, NbHard, pHard, NHard, dHard, tauHard, lambda3D, T_2D_hard, T_3D_hard, fSigma);
	initializeParameters_2(prms_2, kWien, ktWien, NfWien, NsWien, NprWien, NbWien, pWien, NWien, dWien, tauWien, T_2D_wien, T_3D_wien, fSigma);

	//! Declarations
	Video<float> vid, vid_noisy, vid_basic, vid_denoised, vid_diff;
#ifdef OPTICALFLOW
	Video<float> flow;
#endif

	//! Load video
	vid.loadVideo(input_path, firstFrame, lastFrame, frameStep);
	if(kHard == 0)
		vid_basic.loadVideo(inbsc_path, firstFrame, lastFrame, frameStep);

#ifdef OPTICALFLOW
	flow.loadVideo(flow_path, firstFrame, lastFrame-1, frameStep);
#endif

	vid_noisy.resize(vid.sz);
	vid_diff.resize(vid.sz);

	//! Add noise
	if(addnoise)
	{
		cout << endl << "Add noise [sigma = " << fSigma << "] ...";
		VideoUtils::addNoise(vid, vid_noisy, fSigma, verbose);
		cout << "done." << endl;
	}
	else
		vid_noisy = vid;

	//! Denoising
#ifdef OPTICALFLOW
	if (run_vbm3d(fSigma, vid_noisy, flow, vid_basic, vid_denoised, prms_1, prms_2, color_space)
			!= EXIT_SUCCESS)
		return EXIT_FAILURE;
#else
	if (run_vbm3d(fSigma, vid_noisy, vid_basic, vid_denoised, prms_1, prms_2, color_space)
			!= EXIT_SUCCESS)
		return EXIT_FAILURE;
#endif

	//! Compute PSNR and RMSE
	double psnr, rmse;
	double psnr_basic, rmse_basic;
	VideoUtils::computePSNR(vid, vid_basic, psnr_basic, rmse_basic);
	VideoUtils::computePSNR(vid, vid_denoised, psnr, rmse);

	cout << endl << "For noisy video :" << endl;
	cout << "PSNR: " << psnr << endl;
	cout << "RMSE: " << rmse << endl << endl;
	cout << "(basic video) :" << endl;
	cout << "PSNR: " << psnr_basic << endl;
	cout << "RMSE: " << rmse_basic << endl << endl;

	ofstream file(meas_path, ios::out | ios::trunc);
	if(file)
	{
		file << endl << "************" << endl;
		file << "-sigma           = " << fSigma << endl;
		file << "-PSNR_basic      = " << psnr_basic << endl;
		file << "-RMSE_basic      = " << rmse_basic << endl;
		file << "-PSNR            = " << psnr << endl;
		file << "-RMSE            = " << rmse << endl << endl;
		cout << endl;
		file.close();
	}
	else
	{
		cout << "Can't open measures.txt !" << endl;
		return EXIT_FAILURE;
	}

	//! Compute Difference
	cout << endl << "Compute difference...";
	VideoUtils::computeDiff(vid, vid_denoised, vid_diff, fSigma);
	cout << "done." << endl;

	//! save noisy, denoised and differences videos
	cout << endl << "Save videos...";
	vid_noisy.saveVideo(noisy_path, firstFrame, frameStep);
	vid_basic.saveVideo(basic_path, firstFrame, frameStep);
	vid_denoised.saveVideo(final_path, firstFrame, frameStep);
	vid_diff.saveVideo(diff_path, firstFrame, frameStep);

	cout << "done." << endl;

	return EXIT_SUCCESS;
}
