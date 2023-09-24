#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <inttypes.h>

#include "image.h"

int compute_hist(Image *img, uint16_t **dest_hist);
int compute_pix_prob(uint16_t img_size, uint16_t q, uint16_t *hist, float **prob);
int compute_hist_specification(Image *input, Image *spec, Image *out);
int equalize_hist(int q, uint16_t *hist, float *prob);
int equalize_contrast(Image *img);
int image_subsample(Image *img, Image *dst, int factor);
int image_rescale(Image *img, int factor);
int image_requantize(Image *img, int q);
void msb_radixsort(uint16_t *data, int zbin, int obin, uint16_t mask);
void msb_radixsort_index(uint16_t *data, uint16_t *idx, int zbin, int obin, uint16_t mask);
int make_set(uint16_t *data, int n);

#endif // UTIL_H