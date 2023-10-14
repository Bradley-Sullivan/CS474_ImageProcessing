#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>

#include "util.h"

int main(void) {
    Image *input = load_image("catstronaut.pgm");
    
    Image *l = image_laplacian(image_hist_eq(image_gauss(input, 1.4)));

    write_image("laplacian.pgm", l);

    del_image(input); del_image(l);

    return 0;
}






