#ifndef IMAGE_H
#define IMAGE_H

#define IMAGE_DIR "../../images/"

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <ctype.h>
#include <string.h>

typedef struct Image {
    int m, n, q, size;
    uint8_t *data;
} Image;

Image *new_image(int m, int n, int q);
Image *load_image(const char *fname);
int load_header(FILE *fp, Image *img);
int load_data(FILE *fp, Image *img);

int write_image(const char *fname, Image *img);

#endif // IMAGE_H