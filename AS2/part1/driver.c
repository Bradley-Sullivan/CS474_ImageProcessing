#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>
#include <getopt.h>

#include "util.h"

static struct option loptions[] = {
    {"input", required_argument, 0, 'i'},
    {"mask", required_argument, 0, 'm'},
    {"threshold_pix_count", required_argument, 0, 't'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
};

uint32_t calculate_thresh(int pix_lim, uint16_t *histogram, uint16_t q);

int main(int argc, char *argv[]) {
    int option, op_index = 0, pix_lim = 500;
    char in_fname[32], mask_fname[32];
    strcpy(in_fname, "Image.pgm"); strcpy(mask_fname, "Pattern.pgm");

    while ((option = getopt_long(argc, argv, "i:m:t:h", loptions, &op_index)) != -1) {
        switch (option) {
            case 'i':
                if (optarg) strcpy(in_fname, optarg);
                break;
            case 'm':
                if (optarg) strcpy(mask_fname, optarg);
                break;
            case 't':
                if (optarg) sscanf(optarg, "%d", &pix_lim);
                break;
            case 'h':
                printf("Correct Usage:\n\t%s [options]\n", argv[0]);
                printf("\nOPTIONS\n\t-i, --input [string] --> filename of input file (default=Image.pgm)     \
                        \n\t-m, --mask [string] --> filename of mask image (default=Pattern.pgm)            \
                        \n\t-t, --threshold_pixel_count [int] --> max number of pixels to keep in threshold output (default=20)\n");
                return 0;
        }
    }

    Image *input = load_image(in_fname);
    Image *minput = load_image(mask_fname);
    Image *rot_input = rotate_image90(input);
    Image *rot_minput = rotate_image90(minput);
    Mask *mask = mask_image(minput);
    Mask *rot_mask = mask_image(rot_minput);
    
    printf("\nCorrelating %s with mask %s ...\n", in_fname, mask_fname);
    Image *basic_corr = image_correlate(input, mask);

    printf("\nCorrelating rotated %s with rotated mask %s ...\n", in_fname, mask_fname);
    Image *rot_corr = image_correlate(rot_input, rot_mask);

    printf("\nSumming correlations ...\n");
    Image *t = rotate_image90(rot_corr);
    Image *corr_sum = image_add(basic_corr, t);

    write_image("basic_correlation.pgm", basic_corr);
    write_image("rotated_correlation.pgm", rot_corr);
    write_image("correlation_sum.pgm", corr_sum);

    uint16_t *hist = NULL; compute_hist(corr_sum, &hist);

    printf("\nThresholding correlation sum ...\n");
    uint32_t th = calculate_thresh(pix_lim, hist, input->q);
    Image *thresh = image_thresh(corr_sum, th);

    write_image("threshold.pgm", thresh);

    Image *thresh_sum = image_add(thresh, input);
    write_image("thresh_sum.pgm", thresh_sum);

    return 0;
}

uint32_t calculate_thresh(int pix_lim, uint16_t *histogram, uint16_t q) {
    uint32_t pix_ct = 0, i = 0;
    for (i = q - 1; i >= 0; i -= 1) {
        pix_ct += histogram[i];
        if (pix_ct >= pix_lim) break;
    }

    return i;
}

