#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>

#include "util.h"

int main(void) {
    Image *lenna = load_image("Image.pgm");
    Mask *mask = mask_image("Pattern.pgm");

    Image *out = image_correlate(lenna, mask);

    write_image("corr_test.pgm", image_hist_eq(out));

    Image *sum = image_add(lenna, image_hist_eq(out));

    write_image("corr_test_add.pgm", sum);
    write_image("corr_test.pgm", image_hist_eq(out));

    del_image(lenna); del_mask(mask); del_image(out);

    return 0;
}