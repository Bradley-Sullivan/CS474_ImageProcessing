#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "image.h"
#include "util.h"

// #define ARGS
#define NUM_ARGS 1

 /*
    TODO:
        Write functions to compute and operate on image histograms [x]
        Write histogram equalization [x]
        Write histogram specification [x]
        Write image subsampling [x]
        Write image rescaling [x]
        Write image re-quantization [ ]
        Write windowed convolution [ ]
        Write windowed correlation [ ]
        Write local versions of histogram equalizations and specifications
        Write local and global histogram statistics computations
 */

int main(int argc, char *argv[]) {
#ifdef ARGS
    if (argc <= NUM_ARGS) {
        fprintf(stderr, "Error. Incorrect usage.\n  Correct usage: %s [image_filename]\n", argv[0]);
        return 1;
    }
#endif

    // Image *in, *spec, *out;

    // in = load_image("images/boat.pgm");
    // spec = load_image("images/sf.pgm");
    // out = new_image(in->m, in->n, in->q);

    // compute_hist_specification(in, spec, out);

    // write_image("output_4.pgm", out);

    Image *input, *output = NULL;

    input = load_image(argv[1]);

    output = new_image(0,0,0);

    image_subsample(input, output, 8);

    write_image("subsample.pgm", output);

    image_rescale(output, 8);

    write_image("rescale.pgm", output);

    return 0;
}
