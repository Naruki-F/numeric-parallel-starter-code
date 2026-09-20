#define _GNU_SOURCE 
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <sched.h>

#define NUM_THREADS 8
#define RANGE_STEP  100 // n = (idx+1) * 100

typedef struct
{
    int  threadIdx;
    int  n; // total limit  1..n
    long sum;
    long expected; // n(n+1)/2
} threadParams_t;


// POSIX thread declarations and scheduling attributes
//
pthread_t threads[NUM_THREADS];
threadParams_t threadParams[NUM_THREADS];


void *counterThread(void *threadp)
{
    long sum=0;
    int i;
    threadParams_t *p = (threadParams_t *)threadp;

    for(i = 1; i <= p->n; i++)
        sum=sum+i;
    
    p->sum = sum;
    p->expected = ((long)p->n * (p->n + 1)) / 2;
 
    printf("Thread idx=%d on CPU %d: sum[1...%d]=%ld, n(n+1)/2=%ld -> %s\n", p->threadIdx, sched_getcpu(), p->n, p->sum, p->expected, (p->sum == p->expected) ? "PASS" : "FAIL");

    return((void *)0);
}


int main (int argc, char *argv[])
{
   int i;

   for(i=0; i < NUM_THREADS; i++)
   {
       threadParams[i].threadIdx=i;
       threadParams[i].n = (i + 1) * RANGE_STEP; // 100, 200,..., 800

       pthread_create(&threads[i],   // pointer to thread descriptor
                      (void *)0,     // use default attributes
                      counterThread, // thread function entry point
                      (void *)&(threadParams[i]) // parameters to pass in
                     );

   }

   for(i=0;i<NUM_THREADS;i++)
       pthread_join(threads[i], NULL);

   printf("TEST COMPLETE\n");
}
