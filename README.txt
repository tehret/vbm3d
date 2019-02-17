IMPLEMENTATION OF THE VIDEO DENOISING ALGORITHM VBM3D
=====================================================

* Author    : EHRET Thibaud <ehret.thibaud@gmail.com>
* Copyright : (C) 2018 IPOL Image Processing On Line http://www.ipol.im/
* Licence   : GPL v3+, see gpl.txt

OVERVIEW
--------

This source code provides an implementation of VBM3D developped in "Dabov, Kostadin,
Alessandro Foi, and Karen Egiazarian. "Video denoising by sparse 3D transform-domain
collaborative filtering." 2007 15th European Signal Processing Conference. IEEE, 2007".

This code is part of an IPOL (http://www.ipol.im/) publication. Plase cite it
if you use this code as part of your research. (The article is not already published 
at this time)
It is based on Marc Lebrun's code for the image denoising version of this algorithm (BM3D)
available on the BM3D IPOL page (http://www.ipol.im/pub/art/2012/l-bm3d/).

COMPILATION
-----------

The code is compilable on Unix/Linux and hopefully on Mac OS (not tested!). 

Compilation: requires the make program.

Dependencies: FFTW3 and OpenMP (optional). 
For image i/o we use Enric Meinhardt's iio (https://github.com/mnhrdt/iio),
which requires libpng, libtiff and libjpeg.
 
Compile the source code using make.

UNIX/LINUX/MAC:
$ make

Binaries will be created in the current folder.

NOTE: By default, the code is compiled with OpenMP multithreaded
parallelization disabled. If your system supports it you can activate it 
by specifying during the compilation with:
$ make OMP=1
The code will then use the maximum number of thread available on your machine (up to 32). 
This can be reduced by changing the maximum number of threads allowed in lines 126 and 127 
of `vbm3d.cpp`.

USAGE
-----

The following commands have to be run from the current folder:

List all available options:</br>
$ ./VBM3Ddenoising --help

While being a video denoising algorithm, the method takes as input the frames of the video 
and not an actual video. The frames can be extracted using ffmpeg on linux. For example: 
$ ffmpeg -i video.mp4 video/i%04d.png

There is only four mandatory input arguments:
* `-i` the input sequence
* `-f` the index of the first frame
* `-l` the index of the last frame
* `-sigma` the standard deviation of the noise

When providing a sequence that is already noisy the option `-add` should be set to false.

All path should be given using the C standard. For example to reference to the following video:
video/i0000.png
video/i0001.png
video/i0002.png
video/i0003.png
video/i0004.png
video/i0005.png
video/i0006.png
video/i0007.png
video/i0008.png
video/i0009.png

The command for denoising with a noise standard deviation of 20 should be 
$ ./VBM3Ddenoising -i video/i%04d.png -f 0 -l 9 -sigma 20

-----

FILES
-----

This project contains the following source files:
```
	main function:            src/main.cpp
	vbm3d implementation:     src/vbm3d.h
	                          src/vbm3d.cpp
	parameters container:     src/parameters.h
	Basis trans. operations:  src/lib_transforms.h
	                          src/lib_transforms.cpp
	command line parsing:     src/Utilities/cmd_option.h
	image i/o:                src/Utilities/iio.h
	                          src/Utilities/iio.c
	image container:          src/Utilities/LibImages.h
	                          src/Utilities/LibImages.cpp
	video container:          src/Utilities/LibVideoT.hpp
	                          src/Utilities/LibVideoT.cpp
	random number generator:  src/Utilities/mt19937ar.h
	                          src/Utilities/mt19937ar.c
```

ABOUT THIS FILE
---------------

Copyright 2018 IPOL Image Processing On Line http://www.ipol.im/

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.  This file is offered as-is,
without any warranty.
