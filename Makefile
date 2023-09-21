CC = gcc
CFLAGS = -g -Wall -O2

all: pgm

pgm: image.o util.o driver.o
	$(CC) $(CFLAGS) -o pgm image.o util.o driver.o
	rm *.o

driver.o: driver.c image.h
	$(CC) $(CFLAGS) -c driver.c

image.o: image.c image.h
	$(CC) $(CFLAGS) -c image.c

util.o: util.c util.h
	$(CC) $(CFLAGS) -c util.c

clean:
	rm pgm.exe