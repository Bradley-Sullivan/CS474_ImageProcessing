#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <inttypes.h>

#include "image.h"

#define M_PI 3.14159265358979323846

#define PREWITT    1
#define SOBEL      2
#define SCHARR     3

#define ROWMAJOR   1
#define COLMAJOR   -1

#define DEG2RAD    (M_PI / 180.0f)
#define RAD2DEG    (180.0f / M_PI)

int compute_hist(Image *img, uint16_t **dest_hist);
int compute_pix_prob(size_t img_size, uint16_t q, uint16_t *hist, float **prob);
int equalize_hist(int q, uint16_t *hist, float *prob);
float sample_gauss(int x, int y, float sigma);

Image *image_add(Image *a, Image *b);
Image *image_sub(Image *a, Image *b);
Image *image_xor(Image *a, Image *b);
Image *image_and(Image *a, Image *b);

Image *image_thresh(Image *img, uint8_t t);
Image *image_hist_eq(Image *input);
Image *image_specify_hist(Image *input, Image *spec);
Image *image_subsample(Image *input, int factor);
Image *image_rescale(Image *input, int factor);
Image *image_requantize(Image *img, uint8_t bits);
Image *image_correlate(Image *img, Mask *mask);
Image *image_average(Image *img, int k);
Image *image_gauss(Image *img, float sig);
Image *image_median_filter(Image *img, int k);
Image *image_unsharp_mask(Image *img, Image *smooth);
Image *image_high_boost(Image *img, Image *diff, float k);  // need to fix
Image *image_gradient(Image *img, int prewitt_sobel);
Image *image_laplacian(Image *img);

void image_iter_window(Image *data, Image *out, Mask *mask, uint32_t (*op)(uint16_t**, Mask*), int dir);

uint32_t correlate_cb(uint16_t **window, Mask* mask);
uint32_t convolve_cb(uint16_t **window, Mask* mask);
uint32_t laplacian_cb(uint16_t **window, Mask* mask);

void msb_radixsort(uint16_t *data, int zbin, int obin, uint16_t mask);
void msb_radixsort_index(uint16_t *data, uint16_t *idx, int zbin, int obin, uint16_t mask);

uint16_t read_image_window(Image *img, uint16_t *win, uint8_t m, uint8_t n, size_t pos);

#endif // UTIL_H