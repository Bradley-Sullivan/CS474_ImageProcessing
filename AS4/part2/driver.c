#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>
#include <getopt.h>
#include <unistd.h>

#include "util.h"
#include "image.h"

#define SOBEL_X 0
#define SOBEL_Y 1
#define SOBEL_W 3
#define SOBEL_H 3
#define LAPLACIAN_W 3
#define LAPLACIAN_H 3

Image *freq_filter(Image *img, double filter[][4], int fr, int fc);
void visualize_filter(double **d, int m, int n);
double **build_sobel(int dir);
double **build_laplacian();

static struct option loptions[] = {
    {"input", required_argument, 0, 'i'},
    {"sobel", no_argument, 0, 's'},
    {"laplacian", no_argument, 0, 'l'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
};

int main(int argc, char* argv[]) {
    int option, op_index = 0;
    char fname[32]; strcpy(fname, "lenna.pgm");
    int mask = 0;       // default SOBEL
    double sobel_x[SOBEL_H + 1][SOBEL_W + 1] = {
        {0,0,0,0},
        {0,1,0,-1},
        {0,2,0,-2},
        {0,1,0,-1}
    };
    double sobel_y[SOBEL_H + 1][SOBEL_W + 1] = {
        {0,0,0,0},
        {0,1,2,1},
        {0,0,0,0},
        {0,-1,-2,-1}
    };
    double laplacian[LAPLACIAN_H + 1][LAPLACIAN_W + 1] = {
        {0,0,0,0},
        {0,0,1,0},
        {0,1,-4,1},
        {0,0,1,0}
    };

    while ((option = getopt_long(argc, argv, "i:s:l:h", loptions, &op_index)) != -1) {
        switch (option) {
            case 'i':
                if (optarg) strcpy(fname, optarg);
                break;
            case 's':
                mask = 0;
                break;
            case 'l':
                mask = 1;
                break;
            case 'h':
                printf("Correct Usage:\n\t%s [input] [options]\n", argv[0]);
                printf("\nOPTIONS\n\t-i, --input [string] --> filename of input file (default=catstronaut.pgm) \
                        \n\t-s, --sobel     --> filters with x and y Sobel masks                                   \
                        \n\t-l, --laplacian --> filters with Laplacian mask                                    \
                        \n\t-h, --help      --> displays this message");
                return 0;
        }
    }

    Image *input = load_image(fname), *write;

    if (mask) {     // laplacian
        write = freq_filter(input, laplacian, LAPLACIAN_H, LAPLACIAN_W);
        write_image("laplacian.pgm", write); del_image(write);
    } else {        // sobel
        write = freq_filter(input, sobel_x, SOBEL_H, SOBEL_W);
        write_image("sobel_x.pgm", write); del_image(write);
        write = freq_filter(input, sobel_y, SOBEL_H, SOBEL_W);
        write_image("sobel_y.pgm", write); del_image(write);
    }

    del_image(input);

    return 0;
}

Image *freq_filter(Image *img, double filter[][4], int fr, int fc) {
    Image *ret = copy_image(img);
    // creates mask/filter with power of two dimensions for FFT
    Mask *mask = new_mask(is_two_power(img->m) ? img->m : nx_power_two(img->m),   \
                           is_two_power(img->n) ? img->n : nx_power_two(img->n));
    int r = (mask->m / 2) - ((fr + 1) / 2);     // centering row offset
    int c = (mask->n / 2) - ((fc + 1) / 2);     // centering col offset
    double th, **md, **id;
    double p[2];
    uint16_t x, y;

    // insert filter into center of mask
    for (int i = r; i < r + fr + 1; i += 1) {
        for (int k = c; k < c + fc + 1; k += 1) {
            mask->data[i * mask->n + k] = filter[i - r][k - c];
        }
    }

    // obtain complex-data array for mask and image
    md = mask_complex(mask, 1);
    id = image_complex(ret, 1);

    // take mask and image into the frequency domain
    dft2D(id, ret->m, ret->n, -1);
    dft2D(md, mask->m, mask->n, -1);

    // visualize mask's DFT (works for sobel masks, buggy for laplacian)
    visualize_filter(md, mask->m, mask->n);

    // iterate image complex-data array
    for (uint16_t i = 0; i < ret->m; i += 1) {
        for (uint16_t k = 1; k < 2 * ret->n; k += 2) {
            x = ((k - 1) >> 1); y = i;          // compute x,y coord for centering
            // placing the filter into the center of the mask -> multiply by pow(-1,x+y)
            // need to undo it to match up with our image's centered DFT (kinda confusing i know)
            md[i][k] *= pow(-1,x+y);
            md[i][k + 1] *= pow(-1,x+y);
            // apply our mask to our image
            cmult(md[i] + k, id[i] + k, p);
            // store the resultant product
            id[i][k] = p[0];
            id[i][k + 1] = p[1];
        }
    }

    // bring our image back into the spatial domain
    dft2D(id, ret->m, ret->n, 1);

    for (uint16_t i = 0; i < ret->m; i += 1) {
        for (uint16_t k = 1; k < ret->n * 2; k += 2) {
            x = ((k - 1) >> 1); y = i;
            // scales and stores resultant real-component of our image's filtered DFT
            th = CLAMP(0,UINT8_MAX, 1000 * log(1 + (id[i][k] * pow(-1,x+y))));
            ret->data[y * ret->n + x] = th;
        }
    }

    // free our mask
    del_mask(mask);

    return ret;
}

void visualize_filter(double **d, int m, int n) {
    Image *vis = new_image(m, n, UINT8_MAX);
    double im, th;
    uint16_t x, y;
    for (uint16_t i = 0; i < m; i += 1) {  // iterate rows
        for (uint16_t k = 1; k < n * 2; k += 2) {  // ...
            x = ((k - 1) >> 1); y = i;
            im = d[i][k + 1] * pow(-1,x+y);
            // need to be able to visualize negative values (sobel mask)
            // maps complex-component into a "signed" range of UINT8_MAX
            th = 128 - (im * 10000);
            vis->data[i * n + ((k - 1) >> 1)] = CLAMP(0,UINT8_MAX,th);
        }
    }

    write_image("filter.pgm", vis);
    del_image(vis);
}


