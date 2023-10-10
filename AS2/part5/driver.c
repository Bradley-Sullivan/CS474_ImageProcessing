#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>

#include "util.h"

int main(void) {
    Image *input = load_image("catstronaut.pgm");
    Image *input_dx = new_image(input->m, input->n, input->q);
    Image *input_dy = new_image(input->m, input->n, input->q);

    image_gradient(input, input_dx, input_dy, PREWITT);

    Image *psum = image_add(input_dx, input_dy);

    write_image("prewitt_dx.pgm", input_dx);
    write_image("prewitt_dy.pgm", input_dy);
    write_image("prewitt_sum.pgm", psum);

    image_gradient(input, input_dx, input_dy, SOBEL);

    Image *ssum = image_add(input_dx, input_dy);

    write_image("sobel_dx.pgm", input_dx);
    write_image("sobel_dy.pgm", input_dy);
    write_image("sobel_sum.pgm", ssum);

    del_image(input); del_image(input_dx); del_image(input_dy);

    return 0;
}

