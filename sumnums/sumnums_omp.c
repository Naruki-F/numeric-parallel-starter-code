// Exercise #2 1(c): OpenMP version of sumnums.
// 10 threads sum 1..1,000,000 in 10 sub-ranges of 100,000.
// Thread i sums [i*100000+1, (i+1)*100000]; partial sums are added
// after the parallel region and checked against n(n+1)/2.
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#define COUNT       (1000000)
#define NUM_THREADS (10)

long long gsum[NUM_THREADS];

// Current time in seconds (POSIX clock, MONOTONIC_RAW)
static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC_RAW, &t);
    return t.tv_sec + t.tv_nsec / 1e9;
}

int main(int argc, char *argv[])
{
    long range = COUNT / NUM_THREADS; // 100000
    long long gsumall = 0;
    int i, nthreads = 0;

    double t0 = now();

    // Parallel region = thread create + join in one block
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        int  idx   = omp_get_thread_num(); // 0..9
        long start = idx * range + 1;
        long end   = (idx + 1) * range;
        long j;
        long long local = 0;

        for(j = start; j <= end; j++)
            local += j;

        gsum[idx] = local; // own element only

        // sum of [a,b] = b(b+1)/2 - (a-1)a/2
        long long expected = (end * (end + 1)) / 2 - ((start - 1) * start) / 2;
        printf("Thread %d done: sum[%ld...%ld]=%lld, formula=%lld -> %s\n", idx, start, end, local, expected, (local == expected) ? "PASS" : "FAIL");

        #pragma omp single
        nthreads = omp_get_num_threads(); // confirm we got 10
    }

    // Reduction after all threads finished
    for(i = 0; i < NUM_THREADS; i++)
        gsumall += gsum[i];

    long long expected = ((long long)COUNT * (COUNT + 1)) / 2;
    double t1 = now();
    printf("threads=%d gsumall=%lld, n(n+1)/2=%lld -> %s\n", nthreads, gsumall, expected, (gsumall == expected) ? "PASS" : "FAIL");
    printf("start=%.9f end=%.9f elapsed=%.9f s\n", t0, t1, t1 - t0);
    return 0;
}