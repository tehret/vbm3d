//
// Created by Nicola Pierazzo on 20/10/15.
//

#include <cstdlib>
#include <cstring>
#include <fftw3.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include "multiscaler.hpp"

extern "C" {
#include "iio.h"
}

using std::string;

namespace multiscaler {


void dct_inplace(Image &img) {
  int n[] = {img.rows(), img.columns()};
  fftwf_r2r_kind dct2[] = {FFTW_REDFT10, FFTW_REDFT10};
  fftwf_plan plan = fftwf_plan_many_r2r(2,
                                        n,
                                        img.channels(),
                                        img.data(),
                                        NULL,
                                        img.channels(),
                                        1,
                                        img.data(),
                                        NULL,
                                        img.channels(),
                                        1,
                                        dct2,
                                        FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);

#ifdef ISOMETRIC_DCT
  ////> isometric normalization
  // this normalization (and scaling) affects several other functions: 
  //     idct_inplace, decompose, recompose
  // but they can all be removed only by applying the Normalization below 
  double norm_factor = sqrt(.25f / (img.rows() * img.columns()));
  for (int ch = 0; ch < img.channels(); ++ch) {
    for (int row = 0; row < img.rows(); ++row) {
      img.val(0, row, ch) /= sqrt(2.f);
      for (int col = 0; col < img.columns(); ++col) {
        img.val(col, row, ch) *= norm_factor;
      }
    }
    for (int col = 0; col < img.columns(); ++col) {
      img.val(col, 0, ch) /= sqrt(2.f);
    }
  }
#else
  //> Normalization
  for (int i = 0; i < img.rows() * img.columns() * img.channels(); ++i) {
      img.val(i) /= 4 * img.rows() * img.columns();
  }
#endif
}

void idct_inplace(Image &img) {
#ifdef ISOMETRIC_DCT
  ////> isometric normalization
  long double norm_factor = sqrt(.25f / (img.rows() * img.columns()));
  for (int ch = 0; ch < img.channels(); ++ch) {
    for (int row = 0; row < img.rows(); ++row) {
      img.val(0, row, ch) *= sqrt(2.f);
      for (int col = 0; col < img.columns(); ++col) {
        img.val(col, row, ch) *= norm_factor;
      }
    }
    for (int col = 0; col < img.columns(); ++col) {
      img.val(col, 0, ch) *= sqrt(2.f);
    }
  }
#endif

  int n[] = {img.rows(), img.columns()};
  fftwf_r2r_kind idct2[] = {FFTW_REDFT01, FFTW_REDFT01};
  fftwf_plan plan = fftwf_plan_many_r2r(2,
                                        n,
                                        img.channels(),
                                        img.data(),
                                        NULL,
                                        img.channels(),
                                        1,
                                        img.data(),
                                        NULL,
                                        img.channels(),
                                        1,
                                        idct2,
                                        FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
}

const char *pick_option(int *c, char **v, const char *o, const char *d) {
  int id = d ? 1 : 0;
  for (int i = 0; i < *c - id; i++) {
    if (v[i][0] == '-' && 0 == strcmp(v[i] + 1, o)) {
      char *r = v[i + id] + 1 - id;
      for (int j = i; j < *c - id; j++)
        v[j] = v[j + id + 1];
      *c -= id + 1;
      return r;
    }
  }
  return d;
}

Image read_image(const string& filename) {
  int w, h, c;
  float *data = iio_read_image_float_vec(filename.c_str(), &w, &h, &c);
  Image im(data, h, w, c);
  free(data);
  return im;
}

void save_image(const Image& image, const std::string& filename) {
  iio_save_image_float_vec(const_cast<char*>(filename.c_str()),
                           const_cast<float*>(image.data()),
                           image.columns(), image.rows(), image.channels());
}

/*----------------------------------------------------------------------------*/
/* filter an image with a separable kernel. In place.
   Adapted from Rafael Grompone's code
 */
void kernel_convolution(Image& image, float* kernel, int n, float offset)
{
  int X = image.rows();
  int Y = image.columns();

  /* Allocate memory */
  float* tmp = new float[X*Y];

  /* auxiliary variables for the double of the image size */
  int nx2 = 2*X;
  int ny2 = 2*Y;

  for(int c = 0; c < image.channels(); ++c)
  {
      /* x axis convolution */
      for(int x = 0; x < X; ++x)
          for(int y = 0; y < Y; ++y)
          {
              float val = 0.0;
              for(int i = 0; i < n; ++i)
              {
                  int j = x - offset + i;

                  /* symmetry boundary condition */
                  while(j<0) j += nx2;
                  while(j>=nx2) j -= nx2;
                  if( j >= X ) j = nx2-1-j;

                  val += image.val(y,j,c) * kernel[i];
              }
              tmp[x+y*X] = val;
          }

      /* y axis convolution */
      for(int x = 0; x < X; ++x)
          for(int y = 0; y < Y; ++y)
          {
              float val = 0.0;
              for(int i = 0; i < n; ++i)
              {
                  int j = y - offset + i;

                  /* symmetry boundary condition */
                  while(j<0) j += ny2;
                  while(j>=ny2) j -= ny2;
                  if( j >= Y ) j = ny2-1-j;

                  val += tmp[x+j*X] * kernel[i];
              }
              image.val(y,x,c) = val;
          }
  }

  /* free memory */
  delete[] tmp;
}

/* filter an image with a kernel in the given direction.
   Adapted from Rafael Grompone's code
 */
Image directional_kernel_convolution(const Image& image, float* kernel, int n, float offset, int direction)
{
  int X = image.rows();
  int Y = image.columns();

  /* Allocate memory */
  Image out(X, Y, image.channels());

  /* auxiliary variables for the double of the image size */
  int nx2 = 2*X;
  int ny2 = 2*Y;

  for(int c = 0; c < image.channels(); ++c)
  {
      if(direction == 0)
      {
          /* x axis convolution */
          for(int x = 0; x < X; ++x)
              for(int y = 0; y < Y; ++y)
              {
                  float val = 0.0;
                  for(int i = 0; i < n; ++i)
                  {
                      int j = x - offset + i;

                      /* symmetry boundary condition */
                      while(j<0) j += nx2;
                      while(j>=nx2) j -= nx2;
                      if( j >= X ) j = nx2-1-j;

                      val += image.val(y,j,c) * kernel[i];
                  }
                  out.val(y,x,c) = val;
              }
      }
      else
      {
          /* y axis convolution */
          for(int x = 0; x < X; ++x)
              for(int y = 0; y < Y; ++y)
              {
                  float val = 0.0;
                  for(int i = 0; i < n; ++i)
                  {
                      int j = y - offset + i;

                      /* symmetry boundary condition */
                      while(j<0) j += ny2;
                      while(j>=ny2) j -= ny2;
                      if( j >= Y ) j = ny2-1-j;

                      val += image.val(j,x,c) * kernel[i];
                  }
                  out.val(y,x,c) = val;
              }
      }
  }

  return out;
}

/* compute a Gaussian kernel of length n, standard deviation sigma,
   and centered at value mean.

   for example, if mean=0.5, the Gaussian will be centered in the middle point
   between values kernel[0] and kernel[1].

   kernel must be allocated to a size n.

   Adapted from Rafael Grompone's code
 */
void gaussian_kernel(float * kernel, int n, float sigma, float mean)
{
  double sum = 0.;

  /* compute Gaussian kernel */
  for(int i = 0; i < n; ++i)
    {
      int val = (i - mean) / sigma;
      kernel[i] = exp(-0.5*val*val);
      sum += kernel[i];
    }

  /* normalization */
  if(sum > 0.) 
      for(int i = 0; i < n; ++i) 
          kernel[i] /= sum;
}

// Blur the image with a Gaussian kernel of std. dev. sigma
void gblur(Image& image, float sigma)
{
    /* compute the Gaussian kernel */
    /*
       The size of the kernel is selected to guarantee that the first discarded
       term is at least 10^prec times smaller than the central value. For that,
       the half size of the kernel must be larger than x, with
       e^(-x^2/2sigma^2) = 1/10^prec
       Then,
       x = sigma * sqrt( 2 * prec * ln(10) )
       */
    //float prec = 3.0;
    //int offset = (int) ceil( sigma * sqrt( 2.0 * prec * log(10.0) ) );
    //int n = 1 + 2 * offset; /* kernel size */
    //float* kernel = new float[n];
    //gaussian_kernel(kernel, n, sigma, (float) offset);

    //kernel_convolution(image, kernel, n, (float) offset);
    //delete[] kernel;
    float prec = 3.0;
    float offset = ((int) ceil( sigma * sqrt( 2.0 * prec * log(10.0) ) )) -0.5;
    int n = 2 * offset; /* kernel size */
    float* kernel = new float[n];
    gaussian_kernel(kernel, n, sigma, (float) offset);

    kernel_convolution(image, kernel, n, (float) offset);
    delete[] kernel;
}

// Subsample an image by a factor of 2 using a Gaussian kernel
Image gdown(const Image& image)
{
    // Blur the input image using a simple Gaussian convolution
    float sigma = 1.39;
    /* compute the Gaussian kernel */
    /*
       The size of the kernel is selected to guarantee that the first discarded
       term is at least 10^prec times smaller than the central value. For that,
       the half size of the kernel must be larger than x, with
       e^(-x^2/2sigma^2) = 1/10^prec
       Then,
       x = sigma * sqrt( 2 * prec * ln(10) )
       */
    float prec = 3.0;
    int offset = (int) ceil( sigma * sqrt( 2.0 * prec * log(10.0) ) );
    int n = 1 + 2 * offset; /* kernel size */
    float* kernel = new float[n];
    gaussian_kernel(kernel, n, sigma, (float) offset);

    Image blurred = image;
    kernel_convolution(blurred, kernel, n, (float) offset);

    delete[] kernel;
    
    Image output(image.rows()/2, image.columns()/2, image.channels());
    // Subsample the blurred image by a factor of 2
    for (int j = 0; j < output.rows(); ++j)
    for (int k = 0; k < output.columns(); ++k)
    for (int l = 0; l < output.channels(); ++l)
        output.val(k, j, l) = blurred.val(k*2, j*2, l);

    return output;
}

// Upsample an image by a factor of 2 using a Gaussian kernel.
// This function upsamples only by a factor of two. The size is 
// given so to deal with uneven sizes
Image gup(const Image& image, int sx, int sy)
{
    assert(image.rows() == sx/2 && image.columns() == sy/2);
    Image output(sx, sy, image.channels());

    // kernels
    int n = 5; /* kernel size */
    float offset1 = 2;
    float* kernel1 = new float[n];
    gaussian_kernel(kernel1, n, 0.5, offset1);
    float offset2 = 1.5;
    float* kernel2 = new float[n];
    gaussian_kernel(kernel2, n, 0.5, offset2);

    Image temp1 = directional_kernel_convolution(image, kernel1, n, offset1, 0);
    Image temp2 = directional_kernel_convolution(temp1, kernel1, n, offset1, 1);
    // Copy the result at the correct position
    for (int j = 0; j < image.rows(); ++j)
    for (int k = 0; k < image.columns(); ++k)
    for (int l = 0; l < image.channels(); ++l)
        output.val(2*k,2*j,l) = temp2.val(k,j,l);

    temp2 = directional_kernel_convolution(temp1, kernel2, n, offset2, 1);
    // Copy the result at the correct position
    for (int j = 0; j < image.rows(); ++j)
    for (int k = 0; k < image.columns(); ++k)
    for (int l = 0; l < image.channels(); ++l)
        output.val(2*k+1,2*j,l) = temp2.val(k,j,l);

    temp1 = directional_kernel_convolution(image, kernel2, n, offset2, 0);
    temp2 = directional_kernel_convolution(temp1, kernel1, n, offset1, 1);
    // Copy the result at the correct position
    for (int j = 0; j < image.rows(); ++j)
    for (int k = 0; k < image.columns(); ++k)
    for (int l = 0; l < image.channels(); ++l)
        output.val(2*k,2*j+1,l) = temp2.val(k,j,l);

    temp2 = directional_kernel_convolution(temp1, kernel2, n, offset2, 1);
    // Copy the result at the correct position
    for (int j = 0; j < image.rows(); ++j)
    for (int k = 0; k < image.columns(); ++k)
    for (int l = 0; l < image.channels(); ++l)
        output.val(2*k+1,2*j+1,l) = temp2.val(k,j,l);

    // If the requested size are uneven, pad the results
    if(sx % 2 == 1)
    {
        for (int k = 0; k < output.columns(); ++k)
        for (int l = 0; l < output.channels(); ++l)
            output.val(k,sx-1,l) = temp2.val(k,sx-2,l);
    }
    if(sy % 2 == 1)
    {
        for (int j = 0; j < output.columns(); ++j)
        for (int l = 0; l < output.channels(); ++l)
            output.val(sy-1,j,l) = temp2.val(sy-2,j,l);
    }
    if(sx % 2 == 1 && sy % 2 == 1)
    {
        for (int l = 0; l < output.channels(); ++l)
            output.val(sy-1,sx-1,l) = temp2.val(sy-2,sx-1,l);
    }

    delete[] kernel1;
    delete[] kernel2;

    return output;
}

// Lanczos 3 kernel
float lanczos3(float x)
{
    if(std::abs(x) >= 3)
        return 0.;
    if(x == 0.)
        return 1;
    return std::sin(M_PI*x)/(M_PI*x) * std::sin(M_PI*x/3.)/(M_PI*x/3.);
}

void lanczos3_kernel(float * kernel, int n, float offset)
{
  double sum = 0.;

  /* compute Lanczos3 kernel */
  for(int i = 0; i < n; ++i)
    {
      kernel[i] = lanczos3((i-offset)/2.)/2.;
      sum += kernel[i];
    }

  /* normalization */
  if(sum > 0.) 
      for(int i = 0; i < n; ++i) 
          kernel[i] /= sum;
}

// Subsample an image by a factor of 2 using a Lanczos kernel
Image ldown(const Image& image)
{
    /* compute the Lanczos kernel */
    float offset = 6;
    int n = 1 + 2*offset; /* kernel size */
    float* kernel = new float[n];
    lanczos3_kernel(kernel, n, offset);

    Image blurred = image;
    kernel_convolution(blurred, kernel, n, offset);

    delete[] kernel;

    Image output(image.rows()/2, image.columns()/2, image.channels());
    // Subsample the blurred image by a factor of 2
    for (int j = 0; j < output.rows(); ++j)
    for (int k = 0; k < output.columns(); ++k)
    for (int l = 0; l < output.channels(); ++l)
        output.val(k, j, l) = blurred.val(k*2, j*2, l);

    return output;
}

// Upsample an image by a factor of 2 using a Lanczos kernel.
// This function upsamples only by a factor of two. The size is 
// given so to deal with uneven sizes
Image lup(const Image& image, int sx, int sy)
{
    assert(image.rows() == sx/2 && image.columns() == sy/2);
    Image output(sx, sy, image.channels());

    // kernels
    int n = 6; /* kernel size */
    float offset1 = 2.25;
    float* kernel1 = new float[n];
    lanczos3_kernel(kernel1, n, offset1);
    float offset2 = 2.75;
    float* kernel2 = new float[n];
    lanczos3_kernel(kernel2, n, offset2);

    Image temp1 = directional_kernel_convolution(image, kernel1, n, offset1, 0);
    Image temp2 = directional_kernel_convolution(temp1, kernel1, n, offset1, 1);
    // Copy the result at the correct position
    for (int j = 0; j < image.rows(); ++j)
    for (int k = 0; k < image.columns(); ++k)
    for (int l = 0; l < image.channels(); ++l)
        output.val(2*k,2*j,l) = temp2.val(k,j,l);

    temp2 = directional_kernel_convolution(temp1, kernel2, n, offset2, 1);
    // Copy the result at the correct position
    for (int j = 0; j < image.rows(); ++j)
    for (int k = 0; k < image.columns(); ++k)
    for (int l = 0; l < image.channels(); ++l)
        output.val(2*k+1,2*j,l) = temp2.val(k,j,l);

    temp1 = directional_kernel_convolution(image, kernel2, n, offset2, 0);
    temp2 = directional_kernel_convolution(temp1, kernel1, n, offset1, 1);
    // Copy the result at the correct position
    for (int j = 0; j < image.rows(); ++j)
    for (int k = 0; k < image.columns(); ++k)
    for (int l = 0; l < image.channels(); ++l)
        output.val(2*k,2*j+1,l) = temp2.val(k,j,l);

    temp2 = directional_kernel_convolution(temp1, kernel2, n, offset2, 1);
    // Copy the result at the correct position
    for (int j = 0; j < image.rows(); ++j)
    for (int k = 0; k < image.columns(); ++k)
    for (int l = 0; l < image.channels(); ++l)
        output.val(2*k+1,2*j+1,l) = temp2.val(k,j,l);

    // If the requested size are uneven, pad the results
    if(sx % 2 == 1)
    {
        for (int k = 0; k < output.columns(); ++k)
        for (int l = 0; l < output.channels(); ++l)
            output.val(k,sx-1,l) = temp2.val(k,sx-2,l);
    }
    if(sy % 2 == 1)
    {
        for (int j = 0; j < output.columns(); ++j)
        for (int l = 0; l < output.channels(); ++l)
            output.val(sy-1,j,l) = temp2.val(sy-2,j,l);
    }
    if(sx % 2 == 1 && sy % 2 == 1)
    {
        for (int l = 0; l < output.channels(); ++l)
            output.val(sy-1,sx-1,l) = temp2.val(sy-2,sx-1,l);
    }

    delete[] kernel1;
    delete[] kernel2;

    return output;
}

}
