#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>

#include "util.h"
#include "image.h"

#define CLAMP(min,max,val) (val > min ? (val < max ? val : max) : min)

int main(void) {
    Image *img = load_image("lenna.pgm");

    int r, im, x, y;
    double **d = (double**) malloc(sizeof(double*) * img->m);
    for (int i = 0; i < img->m; i += 1) {
        d[i] = (double*) malloc(sizeof(double) * img->n * 2 + 1);
    }

    for (uint16_t i = 0; i < img->n; i += 1) {
        for (uint16_t k = 1; k < img->n * 2 + 1; k += 2) {
            x = ((k - 1) >> 1); y = i;
            d[i][k] = (double)(img->data[y * img->n + x]) * pow(-1, x + y);
            d[i][k + 1] = 0;
        }
    }

    dft2D(d, img->m, img->n, -1);

    for (uint16_t i = 0; i < img->m; i += 1) {
        for (uint16_t k = 1; k < img->n * 2; k += 2) {
            r = d[i][k]; im = d[i][k + 1];
            img->data[i * img->n + ((k - 1) >> 1)] = CLAMP(0, UINT8_MAX, 10 * log(1 + sqrt(r * r + im * im)));
        }
    }

    write_image("test.pgm", img);

    return 0;
}



