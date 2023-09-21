#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "image.h"
#include "util.h"

int main(void) {
    Image *boat = load_image("boat.pgm");
    Image *f16 = load_image("f_16.pgm");
    Image *boat_spec = load_image("sf.pgm");
    Image *f16_spec = load_image("peppers.pgm");
    Image *out = new_image(boat->m, boat->n, boat->q);

    compute_hist_specification(boat, boat_spec, out);
    write_image("boat_spec.pgm", out);
    del_image(out); out = new_image(f16->m, f16->n, f16->q);
    compute_hist_specification(f16, f16_spec, out);
    write_image("f16_spec.pgm", out);
    del_image(out);    
    
    return 0;
}
