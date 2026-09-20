CC=gcc
CFLAGS= -O3 -Wall -fopenmp

all:	primes_seg

primes_seg: primes_seg.c
	$(CC) $(CFLAGS) -o $@ primes_seg.c -lm

clean:
	-rm -f primes_seg *.o
