#include <string>
#include <iostream>
#include <cmath>
#include "multiscaler.hpp"

using namespace multiscaler;
using std::string;
using std::to_string;
using std::cerr;
using std::endl;

int main(int argc, char *argv[]) {
  int ratio = 2.;
  bool usage = pick_option(&argc, argv, "h", nullptr);
  int type = atoi(pick_option(&argc, argv, "t", "0"));
  if (argc != 5 || usage) {
    cerr << "Usage: " << argv[0] << " input prefix levels suffix [-t type: 0=DCT, 1=Gaussian, 2=Lanczos]" << endl;
    exit(EXIT_FAILURE);
  }
  string input = argv[1];
  string output_prefix = argv[2];
  int levels = atoi(argv[3]);
  string output_suffix = argv[4];

  // Read the input image
  Image image = read_image(input);

  // Compute the DCT pyramid
  if(type == 0)
  {
      // DCT of the input image
      dct_inplace(image);

      int h{image.rows()}, w{image.columns()};
      for (int i = 0; i < levels; ++i) {
        // Copy data
        Image output(h, w, image.channels());
#ifdef ISOMETRIC_DCT
        //> isometric normalization DCT scaling (not active in dct_inplace)
        const double scaling = std::sqrt((double)(w*h)/((double)(image.rows()*image.columns())));
#else
        const double scaling = 1.0;
#endif
        for (int j = 0; j < h; ++j) {
          for (int k = 0; k < w; ++k) {
            for (int l = 0; l < image.channels(); ++l) {
              output.val(k, j, l) = image.val(k, j, l) * scaling;
            }
          }
        }

        // Inverse DCT
        idct_inplace(output);

        string filename = output_prefix + to_string(i) + output_suffix;
        save_image(output, filename);

        w /= ratio;
        h /= ratio;
      }
  }
  // Compute the Gaussian pyramid
  else if(type == 1)
  {
      if(ratio != 2.)
      {
          cerr << "Usage: Gaussian pyramid has only been implemented for a ratio (option \"-r\") of 2" << endl;
          exit(EXIT_FAILURE);
      }

      Image current = image;

      string filename = output_prefix + to_string(0) + output_suffix;
      save_image(current, filename);

      for (int i = 1; i < levels; ++i) {
          current = gdown(current);

          string filename = output_prefix + to_string(i) + output_suffix;
          save_image(current, filename);
      }
  }
  // Compute the Lanczos pyramid
  else if(type == 2)
  {
      if(ratio != 2.)
      {
          cerr << "Usage: Lanczos pyramid has only been implemented for a ratio (option \"-r\") of 2" << endl;
          exit(EXIT_FAILURE);
      }
      Image current = image;

      string filename = output_prefix + to_string(0) + output_suffix;
      save_image(current, filename);

      for (int i = 1; i < levels; ++i) {
          current = ldown(current);

          string filename = output_prefix + to_string(i) + output_suffix;
          save_image(current, filename);
      }
  }
  else
  {
    cerr << "Usage: Only possible options available for the type (option \"-t\") of pyramid are DCT (0), Gaussian (1), Lanczos (2)" << endl;
    exit(EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}
