#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "image.h"
#include "util.h"

/* ~experimental~ */
int main(void) {
    Image *lenna = load_image("lenna.pgm");

    image_contrast_eq(lenna);

    Image *lenna_rq_inv = image_requantize(lenna, 32, 1);
    Image *lenna_rq = image_requantize(lenna, 32, 0);

    write_image("lenna_rq_inv.pgm", lenna_rq_inv);
    write_image("lenna_rq.pgm", lenna_rq);
    
    return 0;
}