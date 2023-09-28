#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "image.h"
#include "util.h"

/* ~experimental~ */
int main(void) {
    Image *lenna = load_image("lenna.pgm");

    Image *lenna_rq = image_requantize(lenna, 2, 0);

    write_image("lenna_rq.pgm", lenna_rq);
    
    return 0;
}