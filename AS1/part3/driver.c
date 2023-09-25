#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "image.h"
#include "util.h"

int main(void) {
    Image *boat = load_image("boat.pgm");
    Image *f16 = load_image("f_16.pgm");

    Image *boat_eq = image_contrast_eq(boat);
    Image *f16_eq = image_contrast_eq(f16);

    write_image("boat_eq.pgm", boat_eq);
    write_image("f16_eq.pgm", f16_eq);
    
    return 0;
}
