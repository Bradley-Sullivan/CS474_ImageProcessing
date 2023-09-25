#include "util.h"

// could modify to return pointer to created histogram instead of modifying provided arguments
int compute_hist(Image *img, uint16_t **dest_hist) {
    (*dest_hist) = (uint16_t*) calloc(img->q, sizeof(uint16_t));
    if (!(*dest_hist)) {
        fprintf(stderr, "Failed to allocate memory for histogram.\n");
        return 1;
    }

    for (int i = 0; i < img->size; i += 1) (*dest_hist)[img->data[i]] += 1;

    return 0;
}

// could return pointer to pixel probability map
int compute_pix_prob(size_t img_size, uint16_t q, uint16_t *hist, float **prob) {
    (*prob) = (float*) malloc(sizeof(float) * q);
    if (!(*prob)) {
        fprintf(stderr, "Error allocating memory for pixel probabilities.\n");
        return 1;
    }

    for (int i = 0; i < q; i += 1) {
        (*prob)[i] = (float)hist[i] / img_size;
    }

    return 0;
}

Image *image_specify_hist(Image *input, Image *spec) {
    Image *out = new_image(input->m, input->n, input->q);
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

    return out;
}

int equalize_hist(int q, uint16_t *hist, float *prob) {
    float psum = 0;
    
    for (int i = q; i >= 0; i -= 1) {
        psum = 0;
        for (int k = i; k >= 0; k -= 1) psum += prob[k];
        psum *= q;
        hist[i] = (uint16_t) psum;
    }

    return 0;
}

Image *image_contrast_eq(Image *input) {
    Image *img = copy_image(input);

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

    return img;
}

Image *image_rescale(Image *input, int factor) {
    Image *img = copy_image(input);
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

    return img;
}

Image *image_subsample(Image *input, Image *dst, int factor) {
    if (factor <= 0) {
        fprintf(stderr, "Cannot subsample by nonpositive factor.\n");
        return NULL;
    }

    Image *img = new_image(input->m / factor, input->n / factor, input->q);
    if (!img) return NULL;

    int idx = 0;

    for (int i = 0; i < img->m; i += factor) {
        for (int k = 0; k < img->n; k += factor) {
            img->data[idx++] = img->data[i * img->n + k];
        }
    }    

    return img;
}

Image *image_requantize(Image *input, int q, int inv) {
    Image *img = copy_image(input);
    uint16_t *hist = NULL;

    // compute pixel value frequencies
    compute_hist(img, &hist);

    // deduplicate set of pixel frequencies
    int hsize = make_set(hist, img->q);

    // create array to store indices into histogram
    uint16_t *indices = (uint16_t*) malloc(sizeof(uint16_t) * hsize);
    for (int i = 0; i < hsize; i += 1) indices[i] = i;

    // sort index array on frequency values (preserves pixel->freq. mapping)
    msb_radixsort_index(hist, indices, 0, hsize - 1, 1 << 15);

    for (int i = 0; i < img->size; i += 1) {    // loop image pixels
        uint32_t midx = 0, mdiff = img->q;
        for (int k = fmin(q,hsize) - 1; k >= 0; k -= 1) {   // loop q most frequent pixel values
            if (abs(img->data[i] - indices[k]) < mdiff) {   // compute distance between image pixel v. quantized pixel values
                mdiff = abs(img->data[i] - indices[k]);     // updates new minimum diff
                midx = k;                                   // keeps track of idx of min diff
            }
        }
        img->data[i] = (inv) ? ~midx : midx;                                // updates image pixel with quantized value of least distance
    }

    img->q = q;

    // memory cleanup
    free(hist); free(indices);

    return img;
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

    msb_radixsort_index(data, idx, zbin, oh, mask);     // recursively sorts zero bin
    msb_radixsort_index(data, idx, zh, obin, mask);     // recursively sorts one bin
}

void msb_radixsort(uint16_t *data, int zbin, int obin, uint16_t mask) {
    static uint16_t sw = 0;

    if (mask == 0 || zbin >= obin) return;

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

    msb_radixsort(data, zbin, oh, mask);    // recursively sorts zero bin
    msb_radixsort(data, zh, obin, mask);    // recursively sorts one bin
}

int make_set(uint16_t *data, int n) {
    // NOTE: this takes the first value encountered to be the representative of all duplicate values (frequencies)
    // could change to take median value (median pixel value of duplicate frequencies), swapping/removing others
    int border = n - 1;
    for (int i = 0; i < n; i += 1) {
        for (int k = i; k <= border; k += 1) {
            if (data[i] == data[k]) {
                while (data[k] == data[border]) border -= 1;    // avoids swapping equal values
                int sw = data[k]; data[k] = data[border]; data[border] = sw;
                border -= 1;
            }
        }
    }

    return border;      // return deduplicated partition index
}

Window *new_window(int a, int b, int pos) {
    Window *ret = (Window*) malloc(sizeof(Window));

    ret->w = a;
    ret->a = a; ret->b = b;
    ret->len = a * b;
    ret->pos = pos;

    return ret;
}

Window *new_window_sq(int width, int pos) {
    Window *ret = (Window*) malloc(sizeof(Window));

    ret->w = width;
    ret->a = width; ret->b = width;
    ret->len = width * width;
    ret->pos = pos;

    return ret;
}

// NOTE: wpos specifies index of top-left corner of window in data array
uint16_t *read_window(uint16_t *data, int n, int nwidth, Window *win) {

    // need to take window pos from center coords to top left coords

    // compute data bounds
    int drows = n / nwidth;
    int dcols = nwidth;

    // compute window clipping params
    int tlr = fmax(win->pos / nwidth, 0);
    int brr = fmin(tlr + win->w, drows);
    int tlc = fmax(win->pos % nwidth, 0);
    int brc = fmin(tlc + win->w, dcols);

    // clip window size
    int wrows = brr - tlr;
    int wcols = brc - tlc;

    win->a = wrows;
    win->b = wcols;
    win->len = wrows * wcols;

    uint16_t *wdata = (uint16_t*) malloc(sizeof(uint16_t) * win->len);
    for (int i = 0; i < win->len; i += 1) {
        // compute window-coords
        int r = i / win->b, c = i % win->b;
        // compute window position in data-coords
        int didx = (r + tlr) * nwidth + (tlc + c);
        wdata[i] = data[didx];
    }

    return wdata; 
}




