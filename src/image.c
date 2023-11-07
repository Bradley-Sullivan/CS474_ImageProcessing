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

Image *copy_image(Image* img) {
    Image *ret = new_image(img->m, img->n, img->q);

    if (!memcpy(ret->data, img->data, img->size)) {
        fprintf(stderr, "Error. Failed to copy input data array.\n");
        return NULL;
    }

    return ret;
}

int del_image(Image *img) {
    if (!img) {
        fprintf(stderr, "Cannot free NULL Image pointer.\n");
        return 1;
    } else if (!img->data) {
        fprintf(stderr, "Cannot free NULL Image-data pointer.\n");
        return 1;
    }

    free(img->data);

    free(img);

    return 0;
}

Image *load_image(const char *fname) {
    char fpath[64];    

    sprintf(fpath, "%s%s", IMAGE_DIR, fname);

    Image *img = (Image*) malloc(sizeof(Image));
    FILE *fp = fopen(fpath, "r");

    if (!fp) { 
        fprintf(stderr, "Error opening '%s'.\n", fname); 
        return NULL;
    }

    // gets file size
    fseek(fp, 0, SEEK_END);
    size_t fsize = ftell(fp);
    fseek(fp, 0, SEEK_SET);

    uint8_t *fbuf = (uint8_t*) malloc(sizeof(uint8_t) * fsize);

    // loads entire file into memory
    if (!fread(fbuf, sizeof(uint8_t), fsize, fp)) {
        fprintf(stderr, "Error loading file into buffer.\n");
        return NULL;
    }

    // opens filestream on file buffer
    FILE *fbuf_stream = fmemopen(fbuf, fsize, "r");

    // extracts header information from file buffer
    if (load_header(fbuf_stream, img)) {
        fprintf(stderr, "Error loading file header.\n");
        return NULL;
    }

    if (load_data(fbuf_stream, img)) {
        fprintf(stderr, "Error loading file data.\n");
        return NULL;
    }

    fclose(fp);

    return img;
}

int load_header(FILE *fp, Image *img) {
    char *rbuf = (char*) malloc(sizeof(char) * 64);

    fgets(rbuf, 64, fp);

    if (strcmp(rbuf, "P5\n")) {
        fprintf(stderr, "Error. File is not in PGM format.\n");
        return 1;
    }

    do { fgets(rbuf, 64, fp); } while (rbuf[0] == '#');

    sscanf(rbuf, "%hu %hu\n", &img->m, &img->n);

    fscanf(fp, "%hu\n", &img->q);

    img->size = img->m * img->n;

    return 0;
}

int load_data(FILE *fp, Image *img) {
    uint8_t *img_data = (uint8_t*) malloc(sizeof(uint8_t) * img->size);

    if (!fread(img_data, sizeof(uint8_t), img->size, fp)) {
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
    size_t wsize = fwrite(img->data, sizeof(uint8_t), img->size, fp);
    if (wsize != img->size) {
        fprintf(stderr, "Error writing image data to file.\n");
        return 1;
    }

    fclose(fp);

    return 0;
}

Mask *new_mask(int m, int n) {
    Mask *ret = (Mask*) malloc(sizeof(Mask));

    ret->m = m;
    ret->n = n;
    ret->size = m * n;
    ret->sum = 0;
    ret->data = (double*) malloc(sizeof(double) * ret->size);

    return ret;
}

Mask *copy_mask(Mask *msk) {
    Mask *ret = new_mask(msk->m, msk->n);

    ret->sum = msk->sum;

    if (!memcpy(ret->data, msk->data, msk->size)) {
        fprintf(stderr, "Error. Failed to copy input data array.\n");
        return NULL;
    }

    return ret;
}

Mask *mask_image_file(const char *fname) {
    Image *img = load_image(fname);
    Mask *ret = new_mask(img->m, img->n);

    for (int i = 0; i < img->size; i += 1) {
        ret->sum += img->data[i];
    }

    for (int i = 0; i < img->size; i += 1) {
        ret->data[i] = ((double)img->data[i]) / (double)ret->sum;
    }

    del_image(img);

    return ret;
}

Mask *mask_image(Image *img) {
    Mask *ret = new_mask(img->m, img->n);

    for (int i = 0; i < img->size; i += 1) {
        ret->sum += img->data[i];
    }

    for (int i = 0; i < img->size; i += 1) {
        ret->data[i] = ((double)img->data[i]) / (double)ret->sum;
    }

    return ret;
}

int del_mask(Mask *msk) {
    if (!msk) {
        fprintf(stderr, "Cannot free NULL Mask pointer.\n");
        return 1;
    } else if (!msk->data) {
        fprintf(stderr, "Cannot free NULL Mask-data pointer.\n");
        return 1;
    }

    free(msk->data);

    free(msk);

    return 0;
}
