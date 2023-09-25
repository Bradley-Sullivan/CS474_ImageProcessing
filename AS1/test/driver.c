#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>

#include "util.h"

int main(void) {
    int dwidth = 10;
    int wwidth = 3;

    int dsize = dwidth * dwidth;
    int wsize = wwidth * wwidth;

    uint16_t data[dsize];
    uint16_t *window = NULL;

    srand(time(NULL));
    for (int i = 0; i < dsize; i += 1) {
        data[i] = (uint16_t)rand() % UINT8_MAX;
        if (i % dwidth == 0) printf("\n"); 
        printf(" %d", data[i]);
    }
    printf("\n\n");


    Window *w = new_window_sq(wwidth, 8);
    window = read_window(data, dsize, dwidth, w);

    for (int i = 0; i < w->len; i += 1) {
        printf(" %d", window[i]);
    }
    printf("\n");

    return 0;
}