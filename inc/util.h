#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "image.h"

int compute_hist(Image *img, int **dest_hist);
int compute_pix_prob(int img_size, int q, int *hist, float **prob);
int compute_hist_specification(Image *input, Image *spec, Image *out);
int equalize_hist(int q, int *hist, float *prob);
int equalize_contrast(Image *img);
int image_subsample(Image *img, Image *dst, int factor);
int image_rescale(Image *img, int factor);

#endif // UTIL_H