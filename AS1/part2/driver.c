#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "image.h"
#include "util.h"

int main(void) {
    Image *lenna = load_image("peppers.pgm");

    image_requantize(lenna, 64);

    write_image("lenna_rq.pgm", lenna);
    
    return 0;
}