#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>

#include "util.h"
#include "image.h"

#define CLAMP(min,max,val) (val > min ? (val < max ? val : max) : min)

void recon_phase(Image *img);
void recon_spectrum(Image *img);

int main(void) {
    Image *img = load_image("lenna.pgm");

    recon_phase(img);
    write_image("phase.pgm", img);

    recon_spectrum(img);
    write_image("spectrum.pgm", img);

    del_image(img);
    
    return 0;
}

void recon_phase(Image *img) {
    int x, y;
    double r, im, th;
    double **d = (double**) malloc(sizeof(double*) * img->m);
    for (int i = 0; i < img->m; i += 1) {
        d[i] = (double*) malloc(sizeof(double) * img->n * 2 + 1);
    }

    for (uint16_t i = 0; i < img->m; i += 1) {
        for (uint16_t k = 1; k < img->n * 2 + 1; k += 2) {
            x = ((k - 1) >> 1); y = i;
            d[i][k] = (double)(img->data[y * img->n + x]) * pow(-1, x + y);
            d[i][k + 1] = 0;
        }
    }

    dft2D(d, img->m, img->n, -1);

    for (uint16_t i = 0; i < img->m; i += 1) {
        for (uint16_t k = 1; k < img->n * 2 + 1; k += 2) {
            r = d[i][k]; im = d[i][k + 1];
            th = atan2(im, r);
            d[i][k] = cos(th);
            d[i][k + 1] = sin(th);
        }
    }

    dft2D(d, img->m, img->n, 1);

    for (uint16_t i = 0; i < img->m; i += 1) {
        for (uint16_t k = 1; k < img->n * 2 + 1; k += 2) {
            x = ((k - 1) >> 1); y = i;
            th = CLAMP(0,UINT8_MAX, UINT8_MAX * log(1 + (d[i][k] * pow(-1, x+y))));
            img->data[y * img->n + x] = th;
        }

        free(d[i]);
    }

    free(d);
}

void recon_spectrum(Image *img) {
    int x, y;
    double r, im, th;
    double **d = (double**) malloc(sizeof(double*) * img->m);
    for (int i = 0; i < img->m; i += 1) {
        d[i] = (double*) malloc(sizeof(double) * img->n * 2 + 1);
    }

    for (uint16_t i = 0; i < img->m; i += 1) {
        for (uint16_t k = 1; k < img->n * 2 + 1; k += 2) {
            x = ((k - 1) >> 1); y = i;
            d[i][k] = (double)(img->data[y * img->n + x]) * pow(-1, x + y);
            d[i][k + 1] = 0;
        }
    }

    dft2D(d, img->m, img->n, -1);

    for (uint16_t i = 0; i < img->m; i += 1) {
        for (uint16_t k = 1; k < img->n * 2 + 1; k += 2) {
            r = d[i][k]; im = d[i][k + 1];
            d[i][k] = sqrt(r * r + im * im);
            d[i][k + 1] = 0;
        }
    }

    dft2D(d, img->m, img->n, 1);

    for (uint16_t i = 0; i < img->m; i += 1) {
        for (uint16_t k = 1; k < img->n * 2 + 1; k += 2) {
            x = ((k - 1) >> 1); y = i;
            th = CLAMP(0, UINT8_MAX, d[i][k] * pow(-1, x + y));
            img->data[y * img->n + x] = th;
        }

        free(d[i]);
    }

    free(d);
}


