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
  int type = atoi(pick_option(&argc, argv, "t", "0"));
  if ((argc != 5) || usage) {
    cerr << "Usage: " << argv[0] << " prefix levels suffix output [-c factor] [-t type: 0=DCT, 1=Gaussian, 2=Lanczos]" << endl;
    exit(EXIT_FAILURE);
  }
  string input_prefix = argv[1];
  int levels = atoi(argv[2]);
  string input_suffix = argv[3];
  string output_name = argv[4];

  // Use the bigger image to determine width, height and number of channels
  Image output = read_image(input_prefix + "0" + input_suffix);

  if(type == 0)
  {
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
  }
  // Gaussian
  else if(type == 1)
  {
      Image current = read_image(input_prefix + to_string(levels-1) + input_suffix);
      for (int i = levels-2; i >= 0; --i) {
          // Read level i of the pyramid
          Image image = read_image(input_prefix + to_string(i) + input_suffix);

          Image small = gdown(image);

          // current -= small
          for (int j = 0; j < current.rows(); ++j)
          for (int k = 0; k < current.columns(); ++k)
          for (int l = 0; l < current.channels(); ++l)
              current.val(k,j,l) = current.val(k,j,l) - small.val(k,j,l);

          gblur(current, recompose_factor);

          current = gup(current, image.rows(), image.columns());
          
          // current += image;
          for (int j = 0; j < current.rows(); ++j)
          for (int k = 0; k < current.columns(); ++k)
          for (int l = 0; l < current.channels(); ++l)
              current.val(k,j,l) = current.val(k,j,l) + image.val(k,j,l);
          
          //image + gup(gblur(current - gdown(image),recompose_factor));
      }
      output = current;
  }
  // Lanczos
  else if(type == 2)
  {
      Image current = read_image(input_prefix + to_string(levels-1) + input_suffix);
      for (int i = levels-2; i >= 0; --i) {
          // Read level i of the pyramid
          Image image = read_image(input_prefix + to_string(i) + input_suffix);

          Image small = ldown(image);

          // current -= small
          for (int j = 0; j < current.rows(); ++j)
          for (int k = 0; k < current.columns(); ++k)
          for (int l = 0; l < current.channels(); ++l)
              current.val(k,j,l) = current.val(k,j,l) - small.val(k,j,l);

          gblur(current, recompose_factor);

          current = lup(current, image.rows(), image.columns());
          
          // current += image;
          for (int j = 0; j < current.rows(); ++j)
          for (int k = 0; k < current.columns(); ++k)
          for (int l = 0; l < current.channels(); ++l)
              current.val(k,j,l) = current.val(k,j,l) + image.val(k,j,l);
          
          //image + lup(gblur(current - ldown(image),recompose_factor));
      }
      output = current;//lup(output, 2*output.rows(), 2*output.columns());

  }
  else
  {
    cerr << "Usage: Only possible options available for the type (option \"-t\") of pyramid are DCT (0), Gaussian (1), Lanczos (2)" << endl;
    exit(EXIT_FAILURE);
  }

  // Save the output image
  save_image(output, output_name);

  return EXIT_SUCCESS;
}
