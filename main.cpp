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
,	const int Nf
,	const int Ns
,	const int Npr
,	const int Nb
,	const int p
,	const int N
,	const int d
,	const float tau
,	const float lambda3D
,	const unsigned tau_2D
,	const unsigned tau_3D
,	const float sigma
){
	if(k < 0)	
		prms.k = 8;
	else
		prms.k = k;

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
		prms.d = (7*7)/(255.*255.);
	else
		prms.d = (d*d)/(255.*255.);

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

	if(tau_2D == NONE)
		prms.tau_2D = BIOR;
	else
		prms.tau_2D = tau_2D;

	if(tau_3D == NONE)
		prms.tau_3D = HAAR;
	else
		prms.tau_3D = tau_3D;
}

void initializeParameters_2(
	Parameters& prms
,	const int k
,	const int Nf
,	const int Ns
,	const int Npr
,	const int Nb
,	const int p
,	const int N
,	const int d
,	const float tau
,	const unsigned tau_2D
,	const unsigned tau_3D
,	const float sigma
){
	if(k < 0)	
		prms.k = (sigma > 30) ? 8 : 7;
	else
		prms.k = k;

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
		prms.d = (3*3)/(255.*255.);
	else
		prms.d = (d*d)/(255.*255.);

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

	if(tau_2D == NONE)
		prms.tau_2D = DCT;
	else
		prms.tau_2D = tau_2D;

	if(tau_3D == NONE)
		prms.tau_3D = HAAR;
	else
		prms.tau_3D = tau_3D;
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
	if (argc < 14)
	{
		cout << "usage: BM3D image sigma noisy basic denoised difference bias \
                 difference_bias computeBias tau_2d_hard useSD_hard \
                 tau_2d_wien useSD_wien color_space" << endl;
		return EXIT_FAILURE;
	}

	using std::string;
	const string  input_path = clo_option("-i"    , ""              , "< input sequence");
	const string  inbsc_path = clo_option("-b"    , ""              , "< input basic sequence");
	const string  noisy_path = clo_option("-nisy" , "nisy_%03d.png" , "> noisy sequence");
	const string  final_path = clo_option("-deno" , "deno_%03d.png" , "> denoised sequence");
	const string  basic_path = clo_option("-bsic" , "bsic_%03d.png" , "> basic denoised sequence");
	const string   diff_path = clo_option("-diff" , "diff_%03d.png" , "> difference sequence");
	// TODO: these should be determined automatically from the other outputs.
	const string   bias_path = clo_option("-bdeno", "bdeno_%03d.png", "> bias sequence");
	const string bbasic_path = clo_option("-bbsic", "bbsic_%03d.png", "> bias basic sequence");
	const string  bdiff_path = clo_option("-bdiff", "bdiff_%03d.png", "> bias difference sequence");

	const unsigned firstFrame = clo_option("-f", 0, "first frame");
	const unsigned lastFrame  = clo_option("-l", 0, "last frame");
	const unsigned frameStep  = clo_option("-s", 1, "frame step");

	//! General parameters
	const float fSigma = clo_option("-sigma", 0, "Add noise of standard deviation sigma");
	const bool compute_bias  = (bool) clo_option("-compute-bias", false, "> compute bias outputs");
	const bool verbose  = (bool) clo_option("-verbose"     , true , "> verbose output");
	const unsigned print_prms = (unsigned) clo_option("-print-prms", 0, "> prints parameters for given channels");

	//! VBM3D parameters
	const int kHard = clo_option("-kHard", -1 , "< ");
	const int NfHard = clo_option("-NfHard", -1 , "< ");
	const int NsHard = clo_option("-NsHard", -1 , "< ");
	const int NprHard = clo_option("-NprHard", -1 , "< ");
	const int NbHard = clo_option("-NbHard", -1 , "< ");
	const int pHard = clo_option("-pHard", -1 , "< ");
	const int NHard = clo_option("-NHard", -1 , "< ");
	const int dHard = clo_option("-dHard", -1 , "< ");
	const float tauHard = clo_option("-tauHard", -1. , "< ");
	const int kWien = clo_option("-kWien", -1 , "< ");
	const int NfWien = clo_option("-NfWien", -1 , "< ");
	const int NsWien = clo_option("-NsWien", -1 , "< ");
	const int NprWien = clo_option("-NprWien", -1 , "< ");
	const int NbWien = clo_option("-NbWien", -1 , "< ");
	const int pWien = clo_option("-pWien", -1 , "< ");
	const int NWien = clo_option("-NWien", -1 , "< ");
	const int dWien = clo_option("-dWien", -1 , "< ");
	const float tauWien = clo_option("-tauWien", -1. , "< ");

	const float lambda3D = clo_option("-lambda3d", -1. , "< ");

	//! Variables initialization
	const unsigned tau_2D_hard  = (unsigned) clo_option("-tau2dh", NONE , "< tau_2D_hard");
	if (tau_2D_hard != NONE && tau_2D_hard != DCT && tau_2D_hard != BIOR)
	{
		cout << "tau_2d_hard is not known. Choice is :" << endl;
		cout << " -dct (" << DCT << ")" << endl;
		cout << " -bior (" << BIOR << ")" << endl;
		return EXIT_FAILURE;
	}

	const unsigned tau_2D_wien  = (unsigned) clo_option("-tau2dw", NONE , "< tau_2D_wien");
	if (tau_2D_wien != NONE && tau_2D_wien != DCT && tau_2D_wien != BIOR)
	{
		cout << "tau_2d_wien is not known. Choice is :" << endl;
		cout << " -dct (" << DCT << ")" << endl;
		cout << " -bior (" << BIOR << ")" << endl;
		return EXIT_FAILURE;
	};

	const unsigned tau_3D_hard  = (unsigned) clo_option("-tau3dh", NONE , "< tau_3D_hard");
	if (tau_3D_hard != NONE && tau_3D_hard != DCT && tau_3D_hard != BIOR)
	{
		cout << "tau_3d_hard is not known. Choice is :" << endl;
		cout << " -haar (" << HAAR << ")" << endl;
		cout << " -hadamard (" << HADAMARD << ")" << endl;
		return EXIT_FAILURE;
	}

	const unsigned tau_3D_wien  = (unsigned) clo_option("-tau3dw", NONE , "< tau_3D_wien");
	if (tau_3D_wien != NONE && tau_3D_wien != DCT && tau_3D_wien != BIOR)
	{
		cout << "tau_3d_wien is not known. Choice is :" << endl;
		cout << " -haar (" << HAAR << ")" << endl;
		cout << " -hadamard (" << HADAMARD << ")" << endl;
		return EXIT_FAILURE;
	};

	Parameters prms_1;
	Parameters prms_2;

	initializeParameters_1(prms_1, kHard, NfHard, NsHard, NprHard, NbHard, pHard, NHard, dHard, tauHard, lambda3D, tau_2D_hard, tau_3D_hard, fSigma);
	initializeParameters_2(prms_2, kWien, NfWien, NsWien, NprWien, NbWien, pWien, NWien, dWien, tauWien, tau_2D_wien, tau_3D_wien, fSigma);

	//! Declarations
	Video<float> img, img_noisy, img_basic, img_denoised, img_bias, img_diff;
	Video<float> img_basic_bias;
	Video<float> img_diff_bias;

	//! Load image
	img.loadVideo(input_path, firstFrame, lastFrame, frameStep);
	if(kHard == 0)
		img_basic.loadVideo(inbsc_path, firstFrame, lastFrame, frameStep);

	const unsigned color_space  =  (unsigned) clo_option("-color", 0 , "< color space");

	img_noisy.resize(img.sz);
	img_diff.resize(img.sz);

	//! Add noise
	cout << endl << "Add noise [sigma = " << fSigma << "] ...";
	VideoUtils::addNoise(img, img_noisy, fSigma, verbose);
	cout << "done." << endl;

	//! Denoising
	if (run_vbm3d(fSigma, img_noisy, img_basic, img_denoised, prms_1, prms_2, color_space)
			!= EXIT_SUCCESS)
		return EXIT_FAILURE;

	//! Compute PSNR and RMSE
	double psnr, rmse, psnr_bias, rmse_bias;
	double psnr_basic, rmse_basic, psnr_basic_bias, rmse_basic_bias;
	VideoUtils::computePSNR(img, img_basic, psnr_basic, rmse_basic);
	VideoUtils::computePSNR(img, img_denoised, psnr, rmse);

	cout << endl << "For noisy image :" << endl;
	cout << "PSNR: " << psnr << endl;
	cout << "RMSE: " << rmse << endl << endl;
	cout << "(basic image) :" << endl;
	cout << "PSNR: " << psnr_basic << endl;
	cout << "RMSE: " << rmse_basic << endl << endl;

	ofstream file(bbasic_path, ios::out | ios::trunc);
	if(file)
	{
		file << endl << "************" << endl;
		file << "-sigma           = " << fSigma << endl;
		file << "-PSNR_basic      = " << psnr_basic << endl;
		file << "-RMSE_basic      = " << rmse_basic << endl;
		file << "-PSNR            = " << psnr << endl;
		file << "-RMSE            = " << rmse << endl << endl;
		if (compute_bias)
		{
			file << "-PSNR_basic_bias = " << psnr_basic_bias << endl;
			file << "-RMSE_basic_bias = " << rmse_basic_bias << endl;
			file << "-PSNR_bias       = " << psnr_bias << endl;
			file << "-RMSE_bias       = " << rmse_bias << endl;
		}
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
	VideoUtils::computeDiff(img, img_denoised, img_diff, fSigma);
	cout << "done." << endl;

	//! save noisy, denoised and differences images
	cout << endl << "Save images...";
	img_noisy.saveVideo(noisy_path, firstFrame, frameStep);
	img_basic.saveVideo(basic_path, firstFrame, frameStep);
	img_denoised.saveVideo(final_path, firstFrame, frameStep);
	img_diff.saveVideo(diff_path, firstFrame, frameStep);

	cout << "done." << endl;

	return EXIT_SUCCESS;
}
