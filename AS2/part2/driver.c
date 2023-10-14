#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>
#include <getopt.h>
#include <unistd.h>

#include "util.h"

static struct option loptions[] = {
    {"input", required_argument, 0, 'i'},
    {"average", required_argument, 0, 'a'},
    {"gauss", required_argument, 0, 'g'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
};

int main(int argc, char* argv[]) {
    int option, op_index = 0;
    int avg_size = 3;
    float sigma = 1.4f;
    char fname[32]; strcpy(fname, "catstronaut.pgm");

    while ((option = getopt_long(argc, argv, "i:a:g", loptions, &op_index)) != -1) {
        switch (option) {
            case 'a':
                if (optarg) sscanf(optarg, "%d", &avg_size);
                break;
            case 'g':
                if (optarg) sscanf(optarg, "%f", &sigma);
                break;
            case 'i':
                if (optarg) strcpy(fname, optarg);
                break;
            case 'h':
                printf("Correct Usage:\n\t%s [input] [options]\n", argv[0]);
                printf("\nOPTIONS\n\t-i, --input [string] --> filename of input file (default=catstronaut.pgm)     \
                                \n\t-a, --average [int] --> size of average mask (default=3)\n\t-g, --gauss [float] --> standard deviation (sigma) (default=1.4f)\n\n");
                break;
        }
    }

    Image *input = load_image(fname);
    Image *avg = image_average(input, avg_size);
    Image *gauss = image_gauss(input, sigma);

    write_image("input_avg.pgm", avg);
    write_image("input_gauss.pgm", gauss);

    del_image(input); del_image(avg); del_image(gauss);

    return 0;
}