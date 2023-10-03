#include "util.h"

int compute_hist(Image *img, uint16_t **dest_hist) {
    (*dest_hist) = (uint16_t*) calloc(img->q, sizeof(uint16_t));
    if (!(*dest_hist)) {
        fprintf(stderr, "Failed to allocate memory for histogram.\n");
        return 1;
    }

    for (int i = 0; i < img->size; i += 1) (*dest_hist)[img->data[i]] += 1;

    return 0;
}

int compute_pix_prob(size_t img_size, uint16_t q, uint16_t *hist, float **prob) {
    (*prob) = (float*) malloc(sizeof(float) * q);
    if (!(*prob)) {
        fprintf(stderr, "Error allocating memory for pixel probabilities.\n");
        return 1;
    }

    for (int i = 0; i < q; i += 1) (*prob)[i] = (float)hist[i] / img_size;

    return 0;
}

Image *image_specify_hist(Image *input, Image *spec) {
    Image *out = new_image(input->m, input->n, input->q);
    uint16_t *in_hist = NULL, *spec_hist = NULL;
    float *in_prob = NULL, *spec_prob = NULL;

    // compute input and specified histograms
    compute_hist(input, &in_hist);
    compute_hist(spec, &spec_hist);

    // compute pixel probabilities across input images (PDFs)
    compute_pix_prob(input->size, input->q, in_hist, &in_prob);
    compute_pix_prob(spec->size, spec->q, spec_hist, &spec_prob);

    // obtain equalization transformations (CDFs)
    equalize_hist(input->q, in_hist, in_prob);      // obtain T(r)
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

Image *image_hist_eq(Image *input) {
    Image *img = copy_image(input);

    uint16_t *hist = NULL;
    float *prob = NULL;

    compute_hist(img, &hist);

    // obtain pixel probabilities (PDF)
    compute_pix_prob(img->size, img->q, hist, &prob);

    // compute contrast equalization transform (CDF)
    equalize_hist(img->q, hist, prob);

    // apply equalized transform
    for (int i = 0; i < img->size; i += 1) img->data[i] = hist[img->data[i]];

    // memory cleanup
    free(hist); free(prob);

    return img;
}

Image *image_rescale(Image *input, int factor) {
    int newr = input->m * factor, newc = input->n * factor;
    int newsize = newr * newc;
    
    Image *img = new_image(newr, newc, input->q);

    uint8_t *new_data = (uint8_t*) malloc(sizeof(uint8_t) * newsize);
    for (int i = 0; i < newr; i += 1) {
        for (int k = 0; k < newc; k += 1) {
            int input_idx = (i / factor) * input->n + (k / factor);
            new_data[i * newc + k] = input->data[input_idx];
        }
    }
    
    free(img->data);

    img->data = new_data;

    return img;
}

Image *image_subsample(Image *input, int factor) {
    if (factor <= 0) {
        fprintf(stderr, "Cannot subsample by nonpositive factor.\n");
        return NULL;
    }

    uint16_t m = input->m / factor, n = input->n / factor;
    Image *img = new_image(m, n, input->q);
    if (!img) return NULL;

    int idx = 0;

    for (int i = 0; i < input->m; i += factor) {
        for (int k = 0; k < input->n; k += factor) {
            img->data[idx++] = input->data[i * input->n + k];
        }
    }

    return img;
}

Image *image_requantize(Image *input, int bits, int inv) {
    Image *img = copy_image(input);
    
    uint8_t c = img->q >> (1 << (bits - 1));

    for (int i = 0; i < img->size; i += 1) {
        img->data[i] = (uint8_t)((float)img->data[i] / img->q * (1 << bits) * c);
    }

    return img;
}

// sorts on array indices in-place. preserves index -> value mappings if required
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

// sorts data array in-place
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

uint16_t read_image_window(Image *img, uint16_t *win, uint8_t k, uint16_t pos) {
    int posc = pos % img->n, posr = pos / img->n;

    int tlr = fmax(0, posr - (k >> 1));
    int brr = fmin(img->m, posr + (k >> 1) + 1);
    int tlc = fmax(0, posc - (k >> 1));
    int brc = fmin(img->n, posc + (k >> 1) + 1);

    int wrows = brr - tlr;
    int wcols = brc - tlc;
    int wsize = wrows * wcols;

    int base_offset = tlr * img->n + tlc;
    for (int i = 0; i < wrows; i += 1) {
        for (int k = 0; k < wcols; k += 1) {
            int didx = base_offset + (i * img->n) + k;
            win[i * wcols + k] = img->data[didx];
        }
    }

    return wsize;
}




