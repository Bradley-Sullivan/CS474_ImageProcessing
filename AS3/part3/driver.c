#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>
#include <getopt.h>

#include "util.h"
#include "image.h"

#define CLAMP(min,max,val) (val > min ? (val < max ? val : max) : min)

void recon_phase(Image *img);
void recon_spectrum(Image *img);

static struct option loptions[] = {
    {"input", required_argument, 0, 'i'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
};

int main(int argc, char* argv[]) {
    int option, op_index = 0;
    char fname[32]; strcpy(fname, "catstronaut.pgm");

    while ((option = getopt_long(argc, argv, "i:h", loptions, &op_index)) != -1) {
        switch (option) {
            case 'i':
                if (optarg) strcpy(fname, optarg);
                break;
            case 'h':
                printf("Correct Usage:\n\t%s [input] [options]\n", argv[0]);
                printf("\nOPTIONS\n\t-i, --input [string] --> filename of input file (default=catstronaut.pgm) \
                        \n\t-h, --help --> displays this message.");
                return 0;
        }
    }

    Image *img = load_image(fname);

    // compute image phase-reconstruction
    recon_phase(img);
    write_image("phase.pgm", img);

    // compute image spectrum-reconstruction
    recon_spectrum(img);
    write_image("spectrum.pgm", img);

    del_image(img);
    
    return 0;
}

void recon_phase(Image *img) {
    int x, y;
    double r, im, th;
    // obtain centered complex representation of our data
    double **d = image_complex(img, 1);

    // computes DFT of our complex image data (-1 for forward transform)
    dft2D(d, img->m, img->n, -1);

    for (uint16_t i = 0; i < img->m; i += 1) {
        for (uint16_t k = 1; k < img->n * 2 + 1; k += 2) {
            r = d[i][k]; im = d[i][k + 1];
            // computes phase angle with arctan(im/r) (made simpler with atan2())
            th = atan2(im, r);
            // sets spectrum to 1 while retaining phase information
            d[i][k] = cos(th);
            d[i][k + 1] = sin(th);
        }
    }

    // computes the inverse-DFT of our modified complex-image data (1 for inverse transform)
    dft2D(d, img->m, img->n, 1);

    for (uint16_t i = 0; i < img->m; i += 1) {
        for (uint16_t k = 1; k < img->n * 2 + 1; k += 2) {
            x = ((k - 1) >> 1); y = i;
            // logartihmic "stretching" clamped between 0 and 255
            th = CLAMP(0,UINT8_MAX, 50 * log(1 + (d[i][k] * pow(-1, x+y))));
            img->data[y * img->n + x] = th;
        }

        // memory cleanup
        free(d[i]);
    }

    free(d);
}

void recon_spectrum(Image *img) {
    int x, y;
    double r, im, th;
    // obtain centered complex representation of our data
    double **d = image_complex(img, 1);

    // computes DFT of our complex image data (-1 for forward transform)
    dft2D(d, img->m, img->n, -1);

    for (uint16_t i = 0; i < img->m; i += 1) {
        for (uint16_t k = 1; k < img->n * 2 + 1; k += 2) {
            r = d[i][k]; im = d[i][k + 1];
            // removes phase information by replacing real part with magnitude.
            // phase is calculated with arctan(im/r), if im=0 then phase=0.
            // preserving spectrum information allows us to reconstruct our image
            // with spectrum data only.
            d[i][k] = sqrt(r * r + im * im);
            d[i][k + 1] = 0;
        }
    }

    // computes the inverse-DFT of our modified complex-image data (1 for inverse transform)
    dft2D(d, img->m, img->n, 1);

    for (uint16_t i = 0; i < img->m; i += 1) {
        for (uint16_t k = 1; k < img->n * 2 + 1; k += 2) {
            x = ((k - 1) >> 1); y = i;
            th = CLAMP(0, UINT8_MAX, d[i][k] * pow(-1, x + y));
            img->data[y * img->n + x] = th;
        }

        // memory cleanup
        free(d[i]);
    }

    free(d);
}


