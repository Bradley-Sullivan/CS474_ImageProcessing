#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "image.h"
#include "util.h"

int main(void) {
    Image *boat = load_image("boat.pgm");
    Image *f16 = load_image("f_16.pgm");
    Image *sf = load_image("sf.pgm");
    Image *peppers = load_image("peppers.pgm");

    Image *sf_on_boat = image_specify_hist(boat, sf);

    write_image("sf_on_boat.pgm", sf_on_boat);

    Image *peppers_on_f16 = image_specify_hist(f16, peppers);

    write_image("peppers_on_f16.pgm", peppers_on_f16);
    
    return 0;
}