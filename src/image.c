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
    char fpath[64];    

    sprintf(fpath, "%s%s", IMAGE_DIR, fname);

    Image *img = (Image*) malloc(sizeof(Image));
    FILE *fp = fopen(fpath, "r");

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
    int header_lines = 3;

    for (int i = 0; i < 4; i += 1) {
        fgets(rbuf, 64, fp);
        if (rbuf[0] == '#') { header_lines = 4; break; }
    }

    fseek(fp, 0, SEEK_SET);

    fgets(rbuf, 64, fp);
    if (strcmp(rbuf, "P5\n") != 0) {  // file mode
        fprintf(stderr, "Error. File is not PGM.\n");
        return 1;
    }

    for (int i = 1; i < header_lines; i += 1) {
        fgets(rbuf, 64, fp);
        if (rbuf[0] == '#') {    // comment line
            continue;
        } else {
            img->m = strtol(rbuf, &rhead, 10);
            img->n = atoi(rhead);

            fgets(rbuf, 64, fp);

            img->q = strtol(rbuf, &rhead, 10);

            break;
        }
    }

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
