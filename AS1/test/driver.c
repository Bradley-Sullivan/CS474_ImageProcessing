#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>

#include "util.h"

int main(void) {
    int dwidth = 10;
    int wwidth = 5;

    int dsize = dwidth * dwidth;
    int wsize = wwidth * wwidth;

    uint16_t data[dsize];
    uint16_t *window = NULL;

    srand(time(NULL));
    data[0] = rand() % UINT8_MAX;
    printf(" %d", data[0]);
    for (int i = 1; i < dsize; i += 1) {
        if (i % dwidth == 0) printf("\n");
        data[i] = ~data[i - 1] + 1;   // causes collisions!
        printf(" %d", data[i]);
    }
    printf("\n\n");


    Window *w = new_window(wwidth, 0, wsize);
    window = read_window(data, dsize, dwidth, w);

    for (int i = 0; i < wsize; i += 1) {
        if (!(window[i] + window[(i + wsize) % wsize])) {
            printf(" oob");
        } else {
            printf(" %d", window[i]);
        }
    }
    printf("\n");

    return 0;
}