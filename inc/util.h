#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <inttypes.h>

#include "image.h"

typedef struct Window {
    int w;
    int len;
    int pos;
    int b;
} Window;

int compute_hist(Image *img, uint16_t **dest_hist);
int compute_pix_prob(size_t img_size, uint16_t q, uint16_t *hist, float **prob);
int equalize_hist(int q, uint16_t *hist, float *prob);

Image *image_contrast_eq(Image *input);
Image *image_specify_hist(Image *input, Image *spec);
Image *image_subsample(Image *input, Image *dst, int factor);
Image *image_rescale(Image *input, int factor);
Image *image_requantize(Image *img, int q, int inv);

void msb_radixsort(uint16_t *data, int zbin, int obin, uint16_t mask);
void msb_radixsort_index(uint16_t *data, uint16_t *idx, int zbin, int obin, uint16_t mask);
int make_set(uint16_t *data, int n);

Window *new_window(int width, int pos, int b);
uint16_t *read_window(uint16_t *data, int n, int nwidth, Window *w);
int write_window(uint16_t *data, int n, uint16_t *wdata, Window *w);

#endif // UTIL_H