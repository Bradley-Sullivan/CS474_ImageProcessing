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

float sample_gauss(int x, int y, float sigma) {
    float ex = -( x * x + y * y) / (2 * sigma * sigma);
    return exp(ex) / (2 * M_PI * sigma * sigma);
}

Image *image_add(Image *a, Image *b) {
    Image *result = new_image(a->m, a->n, a->q);
    for (size_t i = 0; i < a->size; i += 1) {
        result->data[i] = (uint8_t) fmin(a->q, a->data[i] + b->data[i]);
    }
    return result;
}

Image *image_sub(Image *a, Image *b) {
    Image *result = new_image(a->m, a->n, a->q);
    for (size_t i = 0; i < a->size; i += 1) {
        result->data[i] = (uint8_t) fmax(0, a->data[i] - b->data[i]);
    }
    return result;
}

Image *image_xor(Image *a, Image *b) {
    Image *result = new_image(a->m, a->n, a->q);
    for (size_t i = 0; i < a->size; i += 1) {
        result->data[i] = (uint8_t) a->data[i] ^ b->data[i];
    }
    return result;
}

Image *image_and(Image *a, Image *b) {
    Image *result = new_image(a->m, a->n, a->q);
    for (size_t i = 0; i < a->size; i += 1) {
        result->data[i] = (uint8_t) a->data[i] & b->data[i];
    }
    return result;
}

Image *image_thresh(Image *img, uint8_t t) {
    Image *ret = copy_image(img);
    for (size_t i = 0; i < img->size; i += 1) {
        ret->data[i] = (ret->data[i] < t) ? 0 : ret->q;
    }
    return ret;
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
            size_t input_idx = (i / factor) * input->n + (k / factor);
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

    size_t idx = 0;
    for (size_t i = 0; i < input->m; i += factor) {
        for (size_t k = 0; k < input->n; k += factor) {
            img->data[idx++] = input->data[i * input->n + k];
        }
    }

    return img;
}

Image *image_requantize(Image *input, uint8_t bits) {
    Image *img = copy_image(input);
    
    int levels = 1 << (bits);           // # of levels
    float comp = (float)input->q / levels;  // dist between levels

    for (int i = 0; i < img->size; i += 1) {
        img->data[i] = (uint8_t) (((float)img->data[i] / UINT8_MAX) * levels) * comp;
    }

    return img;
}

Image *image_correlate(Image *img, Mask *mask) {
    Image *out = new_image(img->m, img->n, img->q);

    image_iter_window(img, out, mask, correlate_cb, ROWMAJOR);

    return out;
}

Image *image_average(Image *img, int k) {
    Image *out = new_image(img->m, img->n, img->q);
    Mask *mask = new_mask(k, k);

    for (int i = 0; i < mask->size; i += 1) {
        mask->data[i] = 1.0f / mask->size;
    }

    image_iter_window(img, out, mask, correlate_cb, ROWMAJOR);

    del_mask(mask);

    return out;
}

Image *image_gauss(Image *img, float sig) {
    Image *out = new_image(img->m, img->n, img->q);
    int k = 5 * sig;
    Mask *kernel = new_mask(k, k);

    // sample gaussian w/ sig
    uint32_t nh = kernel->n >> 1;
    for (int i = 0; i < kernel->n; i += 1) {
        for (int k = 0; k < kernel->n; k += 1) {
            kernel->data[i * kernel->n + k] = sample_gauss(k - nh, i - nh, sig);
        }
    }

    image_iter_window(img, out, kernel, convolve_cb, ROWMAJOR);

    del_mask(kernel);

    return out;
}

Image *image_median_filter(Image *img, int k) {
    int wsize = k * k;
    Image *out = new_image(img->m, img->n, img->q);
    uint16_t *window = (uint16_t*) malloc(sizeof(uint16_t) * wsize);

    for (size_t i = 0; i < img->size; i += 1) {
        read_image_window(img, window, k, k, i);
        msb_radixsort(window, 0, wsize - 1, (uint16_t)0x80);
        out->data[i] = window[wsize >> 1];
    }
    
    return out;
}

Image *image_unsharp_mask(Image *img, Image *smooth) {
    Image *diff = image_sub(img, smooth);
    Image *out = image_add(img, diff);

    del_image(diff);

    return out;
}

Image *image_high_boost(Image *img, Image *diff, float k) {
    Image *out = new_image(img->m, img->n, img->q);

    for (size_t i = 0; i < img->size; i += 1) {
        out->data[i] = img->data[i] + (k * (float)diff->data[i]);
    }

    return out;
}

Image *image_gradient(Image *img, int prewitt_sobel) {
    Image *dx = new_image(img->m, img->n, img->q);
    Image *dy = new_image(img->m, img->n, img->q);
    Image *out = new_image(img->m, img->n, img->q);
    Mask *d = new_mask(3, 3);

    for (int i = 0; i < 3 * 3; i += 1) {
        d->data[i] = 1 - (i % 3);
        if (i / 3 == 1) d->data[i] *= prewitt_sobel;
    }

    image_iter_window(img, dx, d, convolve_cb, ROWMAJOR);

    for (int i = 0; i < 3 * 3; i += 1) {
        d->data[i] = 1 - (i / 3);
        if (i % 3 == 1) d->data[i] *= prewitt_sobel;
    }

    image_iter_window(img, dy, d, convolve_cb, COLMAJOR);

    for (int i = 0; i < img->size; i += 1) {
        out->data[i] = sqrt(dx->data[i] * dx->data[i] + dy->data[i] * dy->data[i]);
    }

    return out;
}

Image *image_laplacian(Image *img) {
    float DD[9] = { 0, 1, 0, 1, -4, 1, 0, 1, 0 };
    Image *out = new_image(img->m, img->n, img->q);
    Mask d = (Mask){3, 3, 9, 0, DD};

    image_iter_window(img, out, &d, convolve_cb, ROWMAJOR);
    
    return out;
}

void image_iter_window(Image *data, Image *out, Mask *mask, uint32_t (*op)(uint16_t**, Mask*), int dir) {
    uint16_t ws_height = (dir >= 0) ? mask->n : mask->m;
    uint16_t ws_width  = (dir >= 0) ? mask->m : mask->n;

    uint16_t **wswap = (uint16_t**) malloc(sizeof(uint16_t*) * ws_height + 1);
    for (int i = 0; i < ws_height + 1; i += 1 ) wswap[i] = (uint16_t*) calloc(ws_width, sizeof(uint16_t));

    uint16_t nh = mask->n >> 1, mh = mask->m >> 1;
    uint16_t i_bound = (dir >= 0) ? data->m : data->n;
    uint16_t k_bound = (dir >= 0) ? data->n : data->m;

    for (uint32_t i = 0; i < i_bound; i += 1) {
        for (uint32_t k = 0; k < k_bound; k += 1) {
            size_t idx = (dir >= 0) ? i * data->n + k : k * data->m + i;
            if (k == 0 && dir >= 0) {
                for (int j = 0; j < nh + 1; j += 1) {
                    read_image_window(data, wswap[j + nh], ws_height, 1, idx + j);
                }
            } else if (k == 0 && dir < 0) {
                for (int j = 0; j < mh + 1; j += 1) {
                    read_image_window(data, wswap[j + mh], 1, ws_width, idx + (j * data->n));
                }
            } else {
                if (dir >= 0) {
                    uint32_t off = nh + (k / (data->n - nh)) * (k - data->n);
                    read_image_window(data, wswap[ws_height], ws_height, 1, idx + off);
                } else {
                    uint32_t off = (mh + (k / (data->m - mh)) * (k - data->m)) * data->n;
                    read_image_window(data, wswap[ws_height], 1, ws_width, idx + off);
                }

                // shift columns ->
                uint16_t *sw = wswap[0];
                for (int j = 0; j < ws_height; j += 1) wswap[j] = wswap[j + 1];
                wswap[ws_height] = sw;
            }

            for (int i = 0; i < mask->n; i += 1) {
                for (int k = 0; k < mask->m; k += 1) {
                    printf(" %hu", wswap[k][i]);
                }
                printf("\n");
            }

            char ch;
            scanf("%c", &ch);

            uint32_t c = op(wswap, mask);
            out->data[idx] = c;
        }
    }

    if (wswap) free(wswap);
}

uint32_t correlate_cb(uint16_t **window, Mask* mask) {
    float acc = 0;
    for (int i = 0; i < mask->size; i += 1) {
        acc += mask->data[i] * (float)window[i % mask->n][i / mask->n];
    }

    return (uint32_t)acc;
}

uint32_t convolve_cb(uint16_t **window, Mask* mask) {
    float acc = 0;
    for (int i = 0; i < mask->size; i += 1) {
        acc += mask->data[i] * (float)window[mask->n - (i % mask->n) - 1][(mask->m) - (i / mask->n) - 1];
    }

    acc = (acc < 0) ? 0 : (acc > UINT32_MAX) ? UINT32_MAX : acc;

    return (uint32_t)acc;
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

    int zh = zbin, oh = obin;

    do {
        if (data[zh] & mask) {
            sw = data[zh]; data[zh] = data[oh]; data[oh] = sw;
            oh -= 1;
        } else {
            zh += 1;
        }
    } while (zh <= oh);

    mask >>= 1;

    msb_radixsort(data, zbin, oh, mask);    // recursively sorts zero bin
    msb_radixsort(data, zh, obin, mask);    // recursively sorts one bin
}

uint16_t read_image_window(Image *img, uint16_t *win, uint8_t m, uint8_t n, size_t pos) {
    uint8_t mh = m >> 1, nh = n >> 1;
    int32_t posc = pos % img->n, posr = pos / img->n;

    int tlr = (posr - mh > 0) ? posr - mh : 0;
    int brr = (posr + mh + 1 < img->m) ? posr + mh + 1 : img->m;
    int tlc = (posc - nh > 0) ? posc - nh : 0;
    int brc = (posc + nh + 1 < img->n) ? posc + nh + 1 : img->n;

    int rowsh = (posr < mh) ? abs(posr - mh) : 0;
    int colsh = (posc < nh) ? abs(posc - nh) : 0;

    int wrows = brr - tlr;
    int wcols = brc - tlc;
    int wsize = wrows * wcols;

    int base_offset = tlr * img->n + tlc;
    for (int i = 0; i < wrows; i += 1) {
        for (int k = 0; k < wcols; k += 1) {
            uint32_t didx = base_offset + (i * img->n) + k;
            win[(i + rowsh) * n + k + colsh] = img->data[didx];
        }
    }

    return wsize;
}






















