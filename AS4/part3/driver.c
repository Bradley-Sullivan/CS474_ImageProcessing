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

void visualize_dft(Image *img, int centered);
void visualize_filter(const char *fname, double **d, int m, int n);
double samp_high_emph(int u, int v, float gh, float gl, float c, float d);
double samp_bw_high_pass(int u, int v, float w, int n);
double samp_bw_band_reject(int u, int v, float w, float r, int cx, int cy, int n);
double samp_bw_notch_reject(int u, int v, float w, int cx, int cy, int n);

static struct option loptions[] = {
    {"input", required_argument, 0, 'i'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
};

int main(int argc, char* argv[]) {
    int option, op_index = 0, c = 0;
    char fname[32]; strcpy(fname, "catstronaut.pgm");

    while ((option = getopt_long(argc, argv, "i:h", loptions, &op_index)) != -1) {
        switch (option) {
            case 'i':
                if (optarg) strcpy(fname, optarg);
                break;
            case 'h':
                printf("Correct Usage:\n\t%s [input] [options]\n", argv[0]);
                printf("\nOPTIONS\n\t-i, --input [string] --> filename of input file (default=catstronaut.pgm) \
                        \n\t-c, --centered --> centers spectrum \
                        \n\t-h, --help --> displays this message");
                return 0;
        }
    }

    Image *input = load_image(fname);
    double **d = image_complex(input, 0), s[2], p[2];
    int x, y;

    for (size_t i = 0; i < input->m; i += 1) {
        for (size_t k = 1; k < 2 * input->n + 1; k += 2) {
            x = ((k - 1) >> 1); y = i;
            d[i][k] = log((double)input->data[y * input->n + x]) * pow(-1,x+y);
        }
    }

    dft2D(d, input->m, input->n, -1);

    float gh, gl, sl, w;
    gh = 1.2;
    gl = 0.35;
    sl = 0.5;
    w = 4;
    for (size_t i = 0; i < input->m; i += 1) {
        for (size_t k = 1; k < 2 * input->n + 1; k += 2) {
            x = (k >> 1) - (input->n >> 1); y = (input->m >> 1) - i;
            s[0] = samp_high_emph(x, y, gh, gl, sl, w);
            s[1] = 0;
            cmult(d[i] + k, s, p);
            d[i][k] = p[0]; d[i][k + 1] = p[1];
        }
    }

    visualize_filter("filtered_dft.pgm", d, input->m, input->n);

    dft2D(d, input->m, input->n, 1);
    
    for (uint16_t i = 0; i < input->m; i += 1) {
        for (uint16_t k = 1; k < input->n * 2 + 1; k += 2) {
            x = ((k - 1) >> 1); y = i;
            // logartihmic "stretching" clamped between 0 and 255
            input->data[y * input->n + x] = 20 * exp(d[i][k] * pow(-1, x+y));
        }
    }

    write_image("out.pgm", input);

    for (size_t i = 0; i < input->m; i += 1) {
        for (size_t k = 1; k < 2 * input->n + 1; k += 2) {
            x = (k >> 1) - (input->n >> 1); y = (input->m >> 1) - i;
            s[0] = samp_high_emph(x, y, gh, gl, sl, w);
            s[1] = 0;
            d[i][k] = s[0];
            d[i][k + 1] = s[1];
        }
    }

    visualize_filter("pure.pgm", d, input->m, input->n);

    return 0;
}

void visualize_filter(const char *fname, double **d, int m, int n) {
    Image *vis = new_image(m, n, 255);
    double r, im, th;

    for (uint16_t i = 0; i < m; i += 1) {  // iterate rows
        for (uint16_t k = 1; k < n * 2; k += 2) {  // ...
            r = d[i][k]; im = d[i][k + 1];
            // logartihmic "stretching" clamped between 0 and 255
            th = CLAMP(0, UINT8_MAX, 255 * log(1 + sqrt(r * r + im * im)));
            vis->data[i * n + ((k - 1) >> 1)] = th;
        }
    }

    write_image(fname, vis);
    del_image(vis);
}

void visualize_dft(Image *img, int centered) {
    double r, im, th;
    // obtain centered complex representation of our data
    double **d = image_complex(img, centered);

    // computes DFT of our complex image data (-1 for forward transform)
    dft2D(d, img->m, img->n, -1);

    for (uint16_t i = 0; i < img->m; i += 1) {  // iterate rows
        for (uint16_t k = 1; k < img->n * 2; k += 2) {  // ...
            r = d[i][k]; im = d[i][k + 1];
            // logartihmic "stretching" clamped between 0 and 255
            th = CLAMP(0, UINT8_MAX, 255 * log(1 + sqrt(r * r + im * im)));
            img->data[i * img->n + ((k - 1) >> 1)] = th;
        }

        // memory cleanup
        free(d[i]);
    }

    free(d);
}

double samp_high_emph(int u, int v, float gh, float gl, float c, float d) {
    float dd = u * u + v * v;
    float gdiff = gh - gl;
    float m = 1 - exp(-c * (dd / (d * d)));

    return gdiff * m + gl;
}

double samp_bw_high_pass(int u, int v, float w, int n) {
    float d = sqrt(u * u + v * v);
    float denom = 1 + pow(w / d, 2 * n);

    return 1 / denom;
}

double samp_bw_band_reject(int u, int v, float w, float r, int cx, int cy, int n) {
    float dd = (u - cx) * (u - cx) + (v - cy) * (v - cy);
    float denom = 1 + pow(((sqrt(dd) * w) / (dd - (r * r))), 2 * n);

    return 1 / denom;
}

double samp_bw_notch_reject(int u, int v, float w, int cx, int cy, int n) {
    float a = samp_bw_high_pass(u - cx, v - cy, w, n);
    float b = samp_bw_high_pass(u + cx, v + cy, w, n);

    return a * b;
}

