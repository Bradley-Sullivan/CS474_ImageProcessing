#include "image.h"

Image *new_image(int m, int n, int q) {
    Image *ret = (Image*)malloc(sizeof(Image));

    ret->m = m;
    ret->n = n;
    ret->q = q;
    ret->size = m * n;

    ret->data = (uint8_t*) calloc(ret->size, sizeof(uint8_t));

    return ret;
}

Image *load_image(const char *fname) {
    Image *img = (Image*) malloc(sizeof(Image));
    FILE *fp = fopen(fname, "r");
    if (!fp) { 
        fprintf(stderr, "Error opening '%s'.\n", fname); 
        return NULL;
    }

    if (load_header(fp, img)) {
        fprintf(stderr, "Error loading file header.\n");
        return NULL;
    }

    if (load_data(fp, img)) {
        fprintf(stderr, "Error loading file data.\n");
        return NULL;
    }

    fclose(fp);

    return img;
}

int load_header(FILE *fp, Image *img) {
    char rbuf[64], *rhead;

    fgets(rbuf, 2, fp);
    if (rbuf[0] != 80 && rbuf[1] != 53) {
        fprintf(stderr, "Error. File is not PGM.\n");
        return 1;
    }

    do {fscanf(fp, "%s\n", rbuf);} while (rbuf[0] == '#');


    img->m = strtol(rbuf, &rhead, 10);
    img->n = atoi(rhead);

    fgets(rbuf, 64, fp);

    img->q = strtol(rbuf, &rhead, 10);

    fscanf(fp, "%d %d\n%d", &img->m, &img->n, &img->q);

    img->size = img->m * img->n;

    return 0;
}

int load_data(FILE *fp, Image *img) {
    uint8_t *img_data = (uint8_t*) malloc(sizeof(uint8_t) * img->size);

    if (!fread(img_data, img->size, 1, fp)) {
        fprintf(stderr, "Error. Image is not %d by %d.\n", img->m, img->n);
        return 1;
    }

    img->data = img_data;

    return 0;
}

int write_image(const char *fname, Image *img) {
    FILE *fp = fopen(fname, "w");
    if (!fp) { fprintf(stderr, "Error opening destination file '%s'.\n", fname); return 1; }

    // writes header
    fprintf(fp, "P5\n%d %d\n%d\n", img->m, img->n, img->q);

    // writes data
    if (!fwrite(img->data, img->size, 1, fp)) {
        fprintf(stderr, "Error writing image data to file.\n");
        return 1;
    }

    fclose(fp);

    return 0;
}
