#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>
#include <getopt.h>
#include <unistd.h>

#include "util.h"
#include "image.h"

#define CLAMP(min,max,val) (val > min ? (val < max ? val : max) : min)

void visualize_dft(Image *img);

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

    Image *input = load_image(fname);
    visualize_dft(input);
    write_image("dft.pgm", input);

    return 0;
}

void visualize_dft(Image *img) {
    int x, y;
    double r, im, th;
    // obtain centered complex representation of our data
    double **d = image_complex(img, 1);

    // computes DFT of our complex image data (-1 for forward transform)
    dft2D(d, img->m, img->n, -1);

    for (uint16_t i = 0; i < img->m; i += 1) {  // iterate rows
        for (uint16_t k = 1; k < img->n * 2; k += 2) {  // ...
            r = d[i][k]; im = d[i][k + 1];
            // logartihmic "stretching" clamped between 0 and 255
            th = CLAMP(0, UINT8_MAX, 20 * log(1 + sqrt(r * r + im * im)));
            img->data[i * img->n + ((k - 1) >> 1)] = th;
        }

        // memory cleanup
        free(d[i]);
    }

    free(d);
}

