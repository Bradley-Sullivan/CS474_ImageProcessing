#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "image.h"
#include "util.h"

int main(void) {
    Image *boat = load_image("boat.pgm");
    Image *f16 = load_image("f_16.pgm");

    equalize_contrast(boat);
    equalize_contrast(f16);

    write_image("boat_eq.pgm", boat);
    write_image("f16_eq.pgm", f16);
    
    return 0;
}
