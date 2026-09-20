#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>

#define COUNT       (1000000)
#define NUM_THREADS (10)

// Note that often the "digit sum" rather than "sum of the digits" is defined as the sum of the digit in each 10's place, but
// that's not what we want to model here.  E.g., Wikipedia - https://en.wikipedia.org/wiki/Digit_sum
//
// What we want to model is an arithmetic series sum or "sum of numbers in range",  best referred to as the "sum or an arithmetic progression":
// * https://en.wikipedia.org/wiki/Arithmetic_progression
//
// This is what we mean summing numbers in the range 1...n, where we know based on series facts that sum(1...n) = n(n+1)/2
//
// We can now have threads sum sub-ranges of a series as a service and then have them add up  the result after a join so that
// sum(1...n) = sum(1...n/2-1) + sum(n/2...n-1) for example.
//
// This sample code provides a simple example of an arithmetic progression sum (sometimes called sum of the digits for simplicity since
// we know sum(0...9)=9(10)/2=45.
//
// Exercise #2 1(b): 10 threads sum 1..1,000,000 in 10 sub-ranges of 100,000.
// Thread i sums [i*100000+1, (i+1)*100000]; main adds the partial sums
// after join and checks the total against n(n+1)/2 = 500000500000.
//
// It should techically be called a sum of a series of numbers in an arithmetic progression.
//

typedef struct
{
    int threadIdx;
    long start;
    long end;
} threadParams_t;


// POSIX thread declarations and scheduling attributes
pthread_t threads[NUM_THREADS];
threadParams_t threadParams[NUM_THREADS];

// Thread specific globals
long long gsum[NUM_THREADS];

void *sumThread(void *threadp)
{
    threadParams_t *p = (threadParams_t *)threadp;
    long long local = 0;
    long i;


    for(i = p->start; i <= p->end; i++)
        local += i; // not prinf inside loop

    gsum[p->threadIdx] = local; // write own element once

    // sum between [a,b] = b(b+1)/2-(a-1)a/2
    long long expected = (p->end * (p->end + 1)) / 2 - ((p->start - 1) * p->start) / 2;

    printf("Thread %d done: sum[%ld...%ld]=%lld, formula=%lld -> %s\n", p->threadIdx, p->start, p->end, local, expected, (local == expected) ? "PASS" : "FAIL");
    return NULL;
}

int main (int argc, char *argv[])
{
   long range = COUNT / NUM_THREADS; // 100000
   long long gsumall = 0;
   int i;

   // initialize gsum array to zero
   for(i=0; i<NUM_THREADS; i++)
       gsum[i]=0;

   printf("Each thread subrange is %ld\n", range);

   for(i=0; i<NUM_THREADS; i++)
   {
      threadParams[i].threadIdx=i;
      threadParams[i].start = i * range + 1; // 1, 100001, ...
      threadParams[i].end = (i + 1) * range; // 100000, 200000, ...

      if(pthread_create(&threads[i], NULL, sumThread, &threadParams[i]) != 0)
        {
            perror("pthread_create");
            exit(-1);
        }
   }

   for(i=0; i<NUM_THREADS; i++)
     pthread_join(threads[i], NULL);

   // Reduction: add the 10 partial sums after all threads have joined
   for(i=0; i<NUM_THREADS; i++)
       gsumall+=gsum[i];

   long long expected = ((long long)COUNT * (COUNT + 1)) / 2; 
   printf("gsumall=%lld, n(n+1)/2=%lld -> %s\n", gsumall, expected, (gsumall == expected) ? "PASS" : "FAIL");

   return 0;
}
