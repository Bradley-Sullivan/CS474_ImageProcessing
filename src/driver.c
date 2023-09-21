#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "image.h"
#include "util.h"

#define NUM_ARGS 1
#define min(a,b) ((a<b)?a:b)

 /*
    TODO:
        Write functions to compute and operate on image histograms [x]
        Write histogram equalization [x]
        Write histogram specification (?)
        Write local versions of histogram equalizations and specifications
        Write local and global histogram statistics computations
 */

typedef struct Option {
    char short_opt;
    const char *long_opt;
    int num_args;
    void (*handler)(const char *arg);
} Option;

// implement specific handlers (not essential)

int main(int argc, char *argv[]) {
    if (argc <= NUM_ARGS) {
        fprintf(stderr, "Error. Incorrect usage.\n  Correct usage: %s [image_filename]\n", argv[0]);
        return 1;
    }

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
