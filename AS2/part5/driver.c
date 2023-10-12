#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>

#include "util.h"

int main(void) {
    Image *input = load_image("catstronaut.pgm");
    Image *sobel_grad = new_image(input->m, input->n, input->q);
    Image *prewitt_grad = new_image(input->m, input->n, input->q);

    prewitt_grad = image_gradient(input, PREWITT);

    sobel_grad = image_gradient(input, SOBEL);

    write_image("sobel_grad.pgm", sobel_grad);
    
    write_image("prewitt_grad.pgm", prewitt_grad);

    del_image(input); del_image(sobel_grad); del_image(prewitt_grad);

    return 0;
}

