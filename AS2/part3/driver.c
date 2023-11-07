#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>

#include "util.h"

void add_snp_noise(Image *img, double pct);

int main(void) {
    Image *input = load_image("lenna.pgm");
    
    add_snp_noise(input, 0.20f);

    Image *out = image_median_filter(input, 3);

    write_image("noise.pgm", input);
    write_image("filtered.pgm", out);

    del_image(input); del_image(out);

    return 0;
}

void add_snp_noise(Image *img, double pct) {
    size_t snp_ct = pct * img->size;
    while (snp_ct--) {
        size_t idx = (size_t)rand() % img->size;
        if (img->data[idx] > (img->q >> 1)) {
            img->data[idx] = 0;
        } else {
            img->data[idx] = img->q;
        }
    }
}
