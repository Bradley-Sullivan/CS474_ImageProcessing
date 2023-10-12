#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>
#include <getopt.h>

#include "util.h"

#define NUM_ARGS 2

int main(int argc, char *argv[]) {
    if (argc - 1 < NUM_ARGS || argc - 1 > NUM_ARGS) {
        fprintf(stderr, "Correct Usage:\n\t%s [input] [mask]\n\n\t[input] --> filename of input image\n\t[mask] --> filename of mask image\n", argv[0]);
        return 1;
    }

    Image *lenna = load_image(argv[1]);

    Mask *mask = mask_image(argv[2]);

    Image *out = image_correlate(lenna, mask);

    write_image("correlation.pgm", out);

    del_image(lenna); del_mask(mask); del_image(out);

    return 0;
}