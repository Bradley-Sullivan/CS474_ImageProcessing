#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>

#include "util.h"

int main(void) {
    Image *input = load_image("catstronaut.pgm");
    
    Image *l = image_laplacian(image_gauss(image_hist_eq(input), 1.4));

    write_image("laplacian.pgm", l);
    write_image("lthresh.pgm", image_thresh(l, 2));

    del_image(input); del_image(l);

    return 0;
}

