#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>

#include "util.h"

int main(void) {
    Image *input = load_image("tools.pgm");
    Image *smooth = image_gauss(input, 3);
    Image *diff = new_image(input->m, input->n, input->q);

    image_sub(input, smooth, diff);

    Image *unsharp = image_unsharp_mask(input, diff);
    Image *boost = image_high_boost(input, diff, 3);

    write_image("unsharp_mask.pgm", unsharp);
    write_image("high_boost.pgm", boost);

    del_image(input); del_image(unsharp); 
    del_image(smooth); del_image(diff);
    del_image(boost);

    return 0;
}


