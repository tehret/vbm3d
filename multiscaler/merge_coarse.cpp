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
  if ((argc != 4) || usage) {
    cerr << "Usage: " << argv[0] << " image coarse result [-c factor]" << endl;
    exit(EXIT_FAILURE);
  }
  string input_name = argv[1];
  string coarse_name = argv[2];
  string output_name = argv[3];

  // Read the fine image
  Image fine = read_image(input_name);

  // DCT of the fine image
  dct_inplace(fine);

  // Read the low frequencies
  Image coarse = read_image(coarse_name);

  // DCT of the low frequencies
  dct_inplace(coarse);

  // Copy data
  for (int j = 0; j < coarse.rows() * recompose_factor; ++j) {
    for (int k = 0; k < coarse.columns() * recompose_factor; ++k) {
      for (int l = 0; l < coarse.channels(); ++l) {
        fine.val(k, j, l) = coarse.val(k, j, l);
      }
    }
  }

  // IDCT of the output image
  idct_inplace(fine);

  save_image(fine, output_name);

  return EXIT_SUCCESS;
}
