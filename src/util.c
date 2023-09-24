#include "util.h"

int compute_hist(Image *img, uint16_t **dest_hist) {
    (*dest_hist) = (uint16_t*) calloc(img->q, sizeof(uint16_t));
    if (!(*dest_hist)) {
        fprintf(stderr, "Failed to allocate memory for histogram.\n");
        return 1;
    }

    for (int i = 0; i < img->m * img->n; i += 1) (*dest_hist)[img->data[i]] += 1;

    return 0;
}

int compute_pix_prob(uint16_t img_size, uint16_t q, uint16_t *hist, float **prob) {
    (*prob) = (float*) malloc(sizeof(float) * q);
    if (!(*prob)) {
        fprintf(stderr, "Error allocating memory for pixel probabilities.\n");
        return 1;
    }

    for (int i = 0; i < q; i += 1) {
        (*prob)[i] = ((float)hist[i]) / img_size;
    }

    return 0;
}

int compute_hist_specification(Image *input, Image *spec, Image *out) {
    uint16_t *in_hist = NULL, *spec_hist = NULL;
    float *in_prob = NULL, *spec_prob = NULL;

    // compute input and specified histograms
    compute_hist(input, &in_hist);
    compute_hist(spec, &spec_hist);

    // compute pixel probabilities across input images
    compute_pix_prob(input->size, input->q, in_hist, &in_prob);
    compute_pix_prob(spec->size, spec->q, spec_hist, &spec_prob);

    // obtain equalization transformations
    equalize_hist(input->q, in_hist, in_prob); // obtain T(r)
    equalize_hist(spec->q, spec_hist, spec_prob);   // obtain G(z)

    // compute inverse mapping for specified histogram
    int *out_map = (int*) malloc(sizeof(int) * input->q);
    for (int i = 0; i < input->q; i += 1) {
        int s = in_hist[i], k = 0;
        
        while (spec_hist[k] < s && k < spec->q) k += 1;
        
        k = (spec_hist[k] > s) ? k - 1 : k;

        out_map[i] = k;
    }

    // apply inverse mapping on input pixel data, write to output image
    for (int i = 0; i < input->size; i += 1) {
        out->data[i] = out_map[input->data[i]];
    }

    // memory cleanup
    free(out_map);
    free(in_hist); free(in_prob);
    free(spec_hist); free(spec_prob);

    return 0;
}

int equalize_hist(int q, uint16_t *hist, float *prob) {
    float psum = 0;
    for (int i = q; i >= 0; i -= 1) {
        psum = 0;
        for (int k = i; k >= 0; k -= 1) psum += prob[k];
        psum *= q;
        hist[i] = (int) psum;
    }

    return 0;
}

int equalize_contrast(Image *img) {
    uint16_t *hist = NULL;
    float *prob = NULL;

    compute_hist(img, &hist);
    
    compute_pix_prob(img->size, img->q, hist, &prob);

    // compute contrast equalization transform
    equalize_hist(img->q, hist, prob);

    // apply equalized transform
    for (int i = 0; i < img->size; i += 1) img->data[i] = hist[img->data[i]];

    // memory cleanup
    free(hist); free(prob);

    return 0;
}

int image_rescale(Image *img, int factor) {
    int newr = img->m * factor, newc = img->n * factor;
    int newsize = newr * newc;

    uint8_t *new_data = (uint8_t*) malloc(sizeof(uint8_t) * newsize);
    for (int i = 0; i < newr; i += 1) {
        for (int k = 0; k < newc; k += 1) {
            int img_idx = (i / factor) * img->n + (k / factor);
            new_data[i * newc + k] = img->data[img_idx];
        }
    }

    free(img->data);

    img->m = newr; img->n = newc;
    img->size = newsize;
    
    img->data = new_data;

    return 0;
}

int image_subsample(Image *img, Image *dst, int factor) {
    if (factor <= 0) {
        fprintf(stderr, "Cannot subsample by nonpositive factor.\n");
        return 1;
    }

    dst->q = img->q;
    dst->m = img->m / factor; dst->n = img->n / factor;
    dst->size = dst->m * dst->n;

    int idx = 0;
    dst->data = (uint8_t*) realloc(dst->data, sizeof(uint8_t) * dst->size);
    for (int i = 0; i < img->m; i += factor) {
        for (int k = 0; k < img->n; k += factor) {
            dst->data[idx++] = img->data[i * img->n + k];
        }
    }    

    return 0;
}

int image_requantize(Image *img, int q) {
    uint16_t *hist = NULL;

    compute_hist(img, &hist);

    int hsize = make_set(hist, img->q);

    printf("hsize: %d\n", hsize);

    uint16_t *indices = (uint16_t*) malloc(sizeof(uint16_t) * hsize);
    for (int i = 0; i < hsize; i += 1) indices[i] = i;

    msb_radixsort_index(hist, indices, 0, hsize - 1, 1 << 15);

    for (int i = 0; i < img->size; i += 1) {
        // int midx = 0, mdiff = img->q;
        // for (int k = fmin(q,hsize) - 1; k >= 0; k -= 1) {
        //     if (abs(img->data[i] - indices[k]) < mdiff) {
        //         mdiff = abs(img->data[i] - indices[k]);
        //         midx = k;
        //     }
        // }
        // img->data[i] = midx;
    }

    img->q = q;

    return 0;
}

void msb_radixsort_index(uint16_t *data, uint16_t *idx, int zbin, int obin, uint16_t mask) {
    static uint16_t sw = 0;

    if (mask == 0 || zbin >= obin) return;

    uint16_t zh = zbin, oh = obin;

    while (zh <= oh) {
        if (data[idx[zh]] & mask) {
            // sw = data[zh]; data[zh] = data[oh]; data[oh] = sw;
            sw = idx[zh]; idx[zh] = idx[oh]; idx[oh] = sw;
            oh -= 1;
        } else {
            zh += 1;
        }
    }

    mask >>= 1;

    msb_radixsort_index(data, idx, zbin, oh, mask);
    msb_radixsort_index(data, idx, zh, obin, mask);
}

void msb_radixsort(uint16_t *data, int zbin, int obin, uint16_t mask) {
    static uint16_t sw = 0;

    if (mask == 0) return;

    uint8_t zh = zbin, oh = obin;

    while (zh <= oh) {
        if (data[zh] & mask) {
            sw = data[zh]; data[zh] = data[oh]; data[oh] = sw;
            oh -= 1;
        } else {
            zh += 1;
        }
    }

    mask >>= 1;

    msb_radixsort(data, zbin, oh, mask);
    msb_radixsort(data, zh, obin, mask);
}

int make_set(uint16_t *data, int n) {
    // NOTE: this takes the first value encountered to be the representative of all duplicate values
    // could change to take median value, swapping/removing others
    int border = n - 1;
    for (int i = 0; i <= border; i += 1) {
        for (int k = i; k <= border; k += 1) {
            if (data[i] == data[k]) {
                while (data[k] == data[border]) border -= 1;
                int sw = data[k]; data[k] = data[border]; data[border] = sw;
                border -= 1;
            }
        }
    }

    return border;
}