//
// Created by Nicola Pierazzo on 20/10/15.
//

#include <cstdlib>
#include <cstring>
#include <fftw3.h>
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

}
