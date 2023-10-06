#ifndef IMAGE_H
#define IMAGE_H

#define IMAGE_DIR "../../images/"

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <ctype.h>
#include <string.h>

typedef struct Image {
    uint16_t m, n, q;
    size_t size;
    uint8_t *data;
} Image;

typedef struct Mask {
    uint16_t m, n;
    size_t size, sum;
    float *data;
} Mask;

Image *new_image(int m, int n, int q);
Image *copy_image(Image* img);
int del_image(Image *img);
Image *load_image(const char *fname);
int load_header(FILE *fp, Image *img);
int load_data(FILE *fp, Image *img);
int write_image(const char *fname, Image *img);

Mask *new_mask(int m, int n);
Mask *copy_mask(Mask *msk);
Mask *mask_image(const char *fname);
int del_mask(Mask *msk);

#endif // IMAGE_H