#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>

#include "util.h"
#include "image.h"

int main(void) {
    Image *lenna = load_image("lenna.pgm");

    Image *lenna_sub = image_subsample(lenna, 8);

    write_image("lenna_sub.pgm", lenna_sub);

    Image *lenna_rsc = image_rescale(lenna_sub, 8);

    write_image("lenna_rsc.pgm", lenna_rsc);

    return 0;
}