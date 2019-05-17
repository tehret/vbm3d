#include <string>
#include <iostream>
#include "multiscaler.hpp"

using namespace multiscaler;
using std::string;
using std::to_string;
using std::cerr;
using std::endl;

int main(int argc, char *argv[]) {
  float recompose_factor = atof(pick_option(&argc, argv, "c", ".8"));
  bool usage = pick_option(&argc, argv, "h", nullptr);
  if ((argc != 5) || usage) {
    cerr << "Usage: " << argv[0] << " prefix levels suffix output [-c factor]" << endl;
    exit(EXIT_FAILURE);
  }
  string input_prefix = argv[1];
  int levels = atoi(argv[2]);
  string input_suffix = argv[3];
  string output_name = argv[4];

  // Use the bigger image to determine width, height and number of channels
  Image output = read_image(input_prefix + "0" + input_suffix);

  // Perform the DCT
  dct_inplace(output);

  for (int i = 1; i < levels; ++i) {
    // Read level i of the pyramid
    Image image = read_image(input_prefix + to_string(i) + input_suffix);

    // Perform the DCT
    dct_inplace(image);

#ifdef ISOMETRIC_DCT
    //> isometric normalization DCT scaling (not active in dct_inplace)
    const double scaling = std::sqrt((double)(output.rows()*output.columns())/((double)(image.rows()*image.columns())));
#else
    const double scaling = 1.0;
#endif
    // Copy data (selected by recompose_factor)
    for (int j = 0; j < image.rows() * recompose_factor; ++j) {
      for (int k = 0; k < image.columns() * recompose_factor; ++k) {
        for (int l = 0; l < image.channels(); ++l) {
          output.val(k, j, l) = image.val(k, j, l) * scaling;
        }
      }
    }
  }

  // IDCT of the output image
  idct_inplace(output);

  // Save the output image
  save_image(output, output_name);

  return EXIT_SUCCESS;
}
