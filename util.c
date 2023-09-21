#include "util.h"

int compute_hist(Image *img, int **dest_hist) {
    (*dest_hist) = (int*) calloc(img->q, sizeof(int));
    if (!(*dest_hist)) {
        fprintf(stderr, "Failed to allocate memory for histogram.\n");
        return 1;
    }

    for (int i = 0; i < img->m * img->n; i += 1) (*dest_hist)[img->data[i]] += 1;

    return 0;
}

int compute_pix_prob(int img_size, int q, int *hist, float **prob) {
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
    int *in_hist = NULL, *spec_hist = NULL;
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

int equalize_hist(int q, int *hist, float *prob) {
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
    int *hist = NULL;
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