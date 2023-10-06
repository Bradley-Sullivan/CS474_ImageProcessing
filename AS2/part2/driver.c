#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>

#include "util.h"

int main(void) {
    Image *input = load_image("bank.pgm");
    Image *avg = image_average(input, 15);
    Image *gauss = image_gauss(input, 3);

    write_image("input_avg.pgm", avg);
    write_image("input_gauss.pgm", gauss);

    del_image(input); del_image(avg); del_image(gauss);

    return 0;
}