Robust Discontinuity Preserving Optical Flow Methods
----------------------------------------------------

*******
SUMMARY
*******

This is a program for optical flow estimation based in the paper: 

Cite: Monzon, N.; Salgado, A.; Sanchez, J., "Regularization Strategies for Discontinuity-Preserving 
Optical Flow Methods," in Image Processing, IEEE Transactions on , vol.PP, no.12, pp.1-1
doi: 10.1109/TIP.2016.2526903
URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7401084&isnumber=4358840

Abstract: 

In this work, we present an implementation of discontinuity-preserving strategies in TV-L 1
optical flow methods. These are based on exponential functions that mitigate the regularization
at image edges, which usually provide precise flow boundaries. Nevertheless, if the smoothing
is not well controlled, it may produce instabilities in the computed motion fields. We present
an algorithm that allows three regularization strategies: The first uses an exponential function
together with a TV process; the second combines this strategy with a small constant that
ensures a minimum isotropic smoothing; the third is a fully automatic approach that adapts
the diffusion depending on the histogram of the image gradients. The last two alternatives
are aimed at reducing the effect of instabilities. In the experiments, we observe that the pure
exponential function is highly unstable while the other strategies preserve accurate motion
contours for a large range of parameters.



The program is part of an IPOL publication:
http://www.ipol.im/pub/algo/mss_optic_flow/

This program is written by 
Nelson Monzón López <nmonzon@ctim.es> CTIM, Universidad de Las Palmas de Gran Canaria
Agustín Salgado de la Nuez <asalgado@dis.ulpgc.es> CTIM, Universidad de Las Palmas de Gran Canaria
Javier Sánchez Pérez <jsanchez@dis.ulpgc.es> CTIM, Universidad de Las Palmas de Gran Canaria 

Version 1, released on February 19, 2016

This software is distributed under the terms of the BSD license (see file license.txt)


***********
COMPILATION
***********

Required environment: Any unix-like system with a standard compilation
environment (make and C and C++ compilers)

Required libraries: libpng, lipjpeg, libtiff

Compilation instructions: run "make" to produce an executable "main" 


*****
USAGE
*****

The program takes two input images, produce an optical flow as output, and
take some parameters.  The meaning of the parameters is thoroughly discussed on
the accompanying IPOL article.


Run:

	./main I1 I2

or

	./main I1 I2 out_file processors method_type alpha gamma lambda nscales fscale zoom_factor TOL inner_iter outer_iter verbose

where:

	I1: first input image
	I2: second input image
	out_file: name of the output optical flow file
	processors: number of threads
	method_type: Integer that selects the regularization strategy     
	alpha: weight of the smoothing term
	gamma: weight of the gradient constancy term
	lambda: It determines the influence of the exponential function in the regularization.
	nscales: desired number of scales
	fscale: finest scale
	zoom_factor: downsampling factor 
	TOL: stopping criterion threshold for the numerical scheme
	inner_iter: number of inner iterations in the numerical scheme
	outer_iter: number of outer iterations in the numerical scheme
	verbose: 0 or 1, for quiet or verbose behaviour

Parameters can be ommited starting from the end of the list, and they will be
assigned reasonable default values.  Examples:

	./main I1.png I2.png flow.flo
	./main I1.png I2.png flow.flo 1 2 145 15 0.5 100 0 0.75 0.0001 1 38 1

If a parameter is given an invalid value it will take the default value. 
If the output filename is omitted, the flow will be saved on a file named "flow.flo".
