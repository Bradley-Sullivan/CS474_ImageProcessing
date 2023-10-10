#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>

#include "util.h"

int main(void) {
    Image *input = load_image("tools.pgm");
    Image *smooth = image_gauss(input, 3);

    Image *unsharp = image_unsharp_mask(input, smooth);
    Image *boost = image_high_boost(input, smooth, 2);

    write_image("unsharp_mask.pgm", unsharp);
    write_image("high_boost.pgm", boost);

    del_image(input); del_image(unsharp); 
    del_image(smooth); del_image(smooth);
    del_image(boost);

    return 0;
}


