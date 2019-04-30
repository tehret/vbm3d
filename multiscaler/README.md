multiscaler
===========

Decompose and recompose an image into a DCT pyramid.

Every layer of the pyramid contains as many DCT (low) frequencies as its number of pixels.
There is no smoothing, so ringing is to be expected. This code is intended to allow multiscale denoising of images.

Building
--------

To compile, use

    $ mkdir build
    $ cd build
    $ cmake .. [-DCMAKE_CXX_COMPILER=/path/of/c++/compiler -DCMAKE_C_COMPILER=/path/of/c/compiler] [-DCMAKE_BUILD_TYPE=Debug]
    $ make

To rebuild, e.g. when the code is modified, use

    $ cd build
    $ make

Using
-----

To decompose an image `input.png` into 3 levels, respectively `pyramid0.tiff`, `pyramid1.tiff` and `pyramid2.tiff`, use

    $ decompose input.png pyramid 3 .tiff

To recompose it (possibly after working on the single layers) on `output.tiff` use

    $ recompose pyramid 3 .tiff output.tiff

To join two layers `fine.tiff` and `coarse.tiff` together, use

    $ merge_coarse fine.tiff coarse.tiff output.tiff

Use the flag `-c alpha` to select the ratio for conservative recomposing.

Extras
------

The files decompose.py, recompose.py, and merge_coarse.py implement the above functions using the wavelets provided by the pywt package.
These python programs use the piio.py package (included) for reading and writing images in the tiff format.

    $ ./decompose.py input.png pyramid 3 .tiff
    $ ./recompose.py pyramid 3 .tiff output.tiff
