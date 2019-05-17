//
// Created by Nicola Pierazzo on 20/10/15.
//

#ifndef MULTISCALER_MULTISCALER_H
#define MULTISCALER_MULTISCALER_H

#include <cstdlib>
#include <string>
#include "Image.hpp"

#define SMART_PARAMETER_INT(n, v) static int n(void)\
{\
  static int smapa_known_ ## n = 0;\
  static int smapa_value_ ## n = v;\
  if (!smapa_known_ ## n) {\
    int r;\
    char *sv = std::getenv(#n);\
    int y;\
    if (sv)\
      r = sscanf(sv, "%d", &y);\
    if (sv && r == 1)\
      smapa_value_ ## n = y;\
    smapa_known_ ## n = 1;\
  }\
  return smapa_value_ ## n;\
}

namespace multiscaler {

// by default don't use isometric DCT for multiscale dct_inplace, 
// idct_inplace, decompose, and  recompose, because it does 
// unneeded extra computations
//#define ISOMETRIC_DCT

void dct_inplace(Image &img);
void idct_inplace(Image &img);
const char *pick_option(int *c, char **v, const char *o, const char *d);
Image read_image(const std::string& filename);
void save_image(const Image& image, const std::string& filename);

}

#endif //MULTISCALER_MULTISCALER_H
