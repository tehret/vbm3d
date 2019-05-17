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
  float ratio = atof(pick_option(&argc, argv, "r", "2."));
  bool usage = pick_option(&argc, argv, "h", nullptr);
  if (argc != 5 || usage) {
    cerr << "Usage: " << argv[0] << " input prefix levels suffix [-r ratio]" << endl;
    exit(EXIT_FAILURE);
  }
  string input = argv[1];
  string output_prefix = argv[2];
  int levels = atoi(argv[3]);
  string output_suffix = argv[4];

  // Read the input image
  Image image = read_image(input);

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

  return EXIT_SUCCESS;
}
