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

    uint16_t x, y;
    double r, c;
    double **fc = (double**) malloc(sizeof(double*) * img->n);
    double **fr = (double**) malloc(sizeof(double*) * img->n);
    for (int i = 0; i < img->n; i += 1) fc[i] = (double*) malloc(sizeof(double) * img->n * 2 + 1);
    for (int i = 0; i < img->n; i += 1) fr[i] = (double*) malloc(sizeof(double) * img->n * 2 + 1);

    for (uint16_t i = 0; i < img->n; i += 1) {
        for (uint16_t k = 1; k < img->n * 2 + 1; k += 2) {
            x = ((k - 1) >> 1); y = i;
            fr[i][k] = (double)(img->data[i * img->n + ((k - 1) >> 1)]) * pow(-1, x + y);
            fr[i][k + 1] = 0;
        }

        fft(fr[i], img->n, -1);
    }

    for (uint16_t i = 0; i < img->n; i += 1) {
        for (uint16_t k = 1; k < img->n * 2 + 1; k += 2) {
            fc[i][k] = fr[(k - 1) >> 1][(i + 1) << 1] / img->n;
            fc[i][k + 1] = fr[(k - 1) >> 1][((i + 1) << 1) + 1] / img->n;
        }

        fft(fc[i], img->n, -1);
    }
    
    for (uint16_t i = 0; i < img->n; i += 1) {
        for (uint16_t k = 1; k < img->n * 2 + 1; k += 2) {
            r = fc[i][k]; c = fc[i][k + 1];
            img->data[i * img->n + ((k - 1) >> 1)] = CLAMP(0, UINT8_MAX, 10 * log(1 + sqrt(r * r + c * c)));
        }
    }

    write_image("test.pgm", img);

    return 0;
}



