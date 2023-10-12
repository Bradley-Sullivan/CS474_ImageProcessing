#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>

#include "util.h"

int main(void) {
    Image *input = load_image("clouds.pgm");
    
    Image *l = image_hist_eq(image_laplacian(image_gauss(input, 1.4)));

    write_image("laplacian.pgm", l);
    write_image("lthresh.pgm", image_thresh(l, 30));

    del_image(input); del_image(l);

    return 0;
}

