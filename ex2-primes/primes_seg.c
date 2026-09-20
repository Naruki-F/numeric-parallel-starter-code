// Exercise #2 Problem 2: primes in [0, 1e9] and semi-prime factoring
// Method: odd-only segmented Sieve of Eratosthenes, parallelized with OpenMP.
//  - Bit k of the bitmap represents the odd number 2k+1 (even numbers are
//    skipped, halving memory: 5e8 bits = 62.5 MB).
//  - Base primes up to sqrt(MAX) are found first with a small simple sieve.
//  - The bitmap is split into segments of whole 64-bit words. Each OpenMP
//    thread sieves complete segments, so no two threads ever write the same
//    word -> no read-modify-write race and no locks needed.
//  - Primes are counted per segment (popcount), a prefix sum gives each
//    segment its offset in the prime list, then the list is filled in parallel.
//  - A semi-prime SP = P1*P2 is factored by trial division with the prime list
//    (P1 <= sqrt(SP)), then P2 = SP/P1 is checked against the list.
// Usage: ./primes_seg [threads]
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define MAX (1000000000ULL) // upper limit of the range for searching for prime numbers (0 to 1 billion). ULL = unsigned long long
#define NBITS (MAX / 2) // number of bits in a bitmap. Since I only deal with odd numbers, it is half of MAX
#define NWORDS (NBITS / 64) // number of units when expressed in 64-bit units (words). NBITS/64
#define SEG_WORDS (4096) // number of words per segment

static uint64_t *bitmap; // bitmap itself. A 1 indicates a prime number; a 0 indicates a composite number
static uint32_t *primelist; // all primes in ascending order
static uint64_t  nprimes; // number of prime numbers

static double now(void)
{
    struct timespec t; // struct timespec t; = tv_sec and tv_nsec
    clock_gettime(CLOCK_MONOTONIC_RAW, &t); // read the current time from the specified type of clock and writes it to 't'

    return t.tv_sec + t.tv_nsec/1e9; // tv_nesec is nanosecond so it's divided by 10^9
}

// Binary search: x in the prime list or not
static int in_list(uint64_t x)
{
    // lo, hi, and mid are index of list
    int64_t lo = 0, hi = (int64_t)nprimes - 1; 
    while(lo <= hi)
    {
        int64_t mid = (lo + hi)/2; // mid index in the list
        if(primelist[mid] == x) return 1; // find x in the prime list
        if(primelist[mid] < x) lo = mid + 1; else hi = mid - 1; // if p[mid] <= x, the answer is to the right
    }

    return 0;
}

// Index of the first prime in the list that is > x (binary search)
static int64_t upper_index(uint64_t x)
{
    // lo, hi, and mid are index of list
    int64_t lo = 0, hi = (int64_t)nprimes;
    while(lo < hi)
    {
        int64_t mid = (lo + hi)/2; // mid index in the list
        if(primelist[mid] <= x) lo = mid + 1; else hi = mid; // if p[mid] < x means, x is to the right of the mid, otherwise, x is to the left
    }
    return lo;
}

// Find P1 <= P2 with SP = P1*P2 using the prime list. Returns 1 if found.
// If SP = P1 x P2 (where P1 ≤ P2), then P1 must be less than or equal to square root of SP.
static int factor_sp(uint64_t sp, uint64_t *p1, uint64_t *p2)
{
    // integer square root of sp (P1 <= sqrt(SP) for the smaller factor)
    uint64_t limit = (uint64_t)sqrtl((long double)sp);
    // The largest integer that satisfies limit^2 ≤ SP should be adjusted to match limit
    while(limit * limit > sp) limit--;
    while((limit + 1) * (limit + 1) <= sp) limit++;

    int64_t ub = upper_index(limit); // determine the scope of the investigation: only test primes <= sqrt(SP)
    int64_t found = INT64_MAX; // smallest index of a prime divisor. INT64_MAX = still not find it
    int64_t i;

    // parallel for: loop btw 0 and ub is distributed among threads
    // schedule(static): Divides the range into equal, contiguous blocks based on the number of threads
        // Ex. 4 threads, each block covers 1/4 of the range from the beginning
    // reduction(min:found): Each thread maintains its own copy of found. The initial value is the identity element of `min`, which is `INT64_MAX` 
        // Since each thread modifies only its own copy, there are no concurrent writes to shared variables (data races)
        // When the loop ends, OpenMP sets the original `found` to the minimum value among all threads' copies.
    #pragma omp parallel for schedule(static) reduction(min:found)
    for(i = 0; i < ub; i++)
        if(sp % primelist[i] == 0 && i < found) // if a prime number that can be divided is found && if it is in the smallest position
            found = i;

    if(found == INT64_MAX) return 0; // if not find a prime number

    *p1 = primelist[found];
    *p2 = sp / *p1;

    // prime check: P2 must also be a prime < MAX and in the list
    return (*p2 < MAX) && in_list(*p2);
}

int main(int argc, char *argv[])
{
    // Setting the Number of Threads
    int nthreads = (argc > 1) ? atoi(argv[1]) : omp_get_max_threads(); // if ./primes_seg 4, number (4) will be set as the number of threads
    if(nthreads < 1) { fprintf(stderr, "threads must be >= 1\n"); exit(-1); }
    omp_set_num_threads(nthreads); // omp_set_num_threads to specify the number of threads for subsequent parallel processing

    double t_start = now();

    // Base primes up to sqrt(MAX) with a small simple sieve
    uint32_t root = (uint32_t)sqrt((double)MAX) + 1; // sqrt 1 billion = 31,622.77… When converted to an integer, this becomes 31,622, and adding 1 gives 31,623
    // calloc(numberOfElements, sizePerElement), and all elements are initialized to 0. 0 = prime candidate, 1 = composite number
    char *small = calloc(root + 1, 1); // array that uses one byte for each number from 0 to 31,623
    uint32_t *base = malloc(sizeof(uint32_t) * root); // array that stores the prime numbers found with length of root
    uint32_t nbase = 0, a, b; // number of prime numbers found
    if(!small || !base) { perror("malloc base"); exit(-1); }
    
    // look only at the odd numbers: 3, 5, 7, 9, ...
    for(a = 3; a <= root; a += 2)
        if(!small[a]) // if not yet crossed out → prime number
        {
            base[nbase++] = a; // odd base primes only
            for(b = a * a; b <= root; b += 2 * a) small[b] = 1; // remove multiples of a
        }

    // Bitmap: all odd numbers start as "prime"
    bitmap = malloc(NWORDS * sizeof(uint64_t)); // 7,812,500 words x 8 bytes = 62.5 MB
    if(!bitmap) { perror("malloc bitmap"); exit(-1); }

    uint64_t nseg = (NWORDS + SEG_WORDS - 1) / SEG_WORDS; // 7,812,500/4,096 = 1,907.3... so the number of intervals is 1,908
    uint64_t *segcnt = calloc(nseg + 1, sizeof(uint64_t)); // array to store the number of prime numbers in each interval, initialized to 0
    if(!segcnt) { perror("malloc segcnt"); exit(-1); }
    int64_t s;

    double t_sieve = now();
    // 1,908 segments s are distributed among threads
    // schedule(dynamic, 4): take 4 intervals at a time, in the order that threads become available
    #pragma omp parallel for schedule(dynamic, 4)
    for(s = 0; s < (int64_t)nseg; s++)
    {
        // set the range of the interval and initialize
        // interval s covers the range from word w0 to, but not including, word w1
        uint64_t w0 = s * SEG_WORDS;
        uint64_t w1 = (w0 + SEG_WORDS < NWORDS) ? w0 + SEG_WORDS : NWORDS; // In the final segment, since w0 + 4096 exceeds NWORDS, so stop at NWORDS
        memset(&bitmap[w0], 0xFF, (w1 - w0) * sizeof(uint64_t)); // use `memset` to set all bytes in the range to 0xFF (all bits set to 1), and start under the assumption that “all odd numbers in the range are prime.”

        // range of numbers represented by the interval
        uint64_t lo = w0 * 64 * 2 + 1; // first odd number in segment
        uint64_t hi = w1 * 64 * 2 + 1; // one past last odd number
        uint32_t k;
        
        // eliminate multiples of prime numbers 
        // use the prime numbers in ascending order; 
            // once p^2 is greater than or equal to hi, there are no more primes to remove in this interval, so exit the loop
        for(k = 0; k < nbase; k++)
        {
            uint64_t p = base[k];
            if(p * p >= hi) break;
            uint64_t m = ((lo + p - 1) / p) * p; // first multiple >= lo
            if(m < p * p) m = p * p; // smaller ones done by smaller primes
            if((m & 1) == 0) m += p; // odd multiples only. if even number, add p to become odd num
            
            // move forward in increments of 2p, deleting only the multiples of odd numbers. If p = 3, then 524,289 -> 524,295 -> 524,301
            for(; m < hi; m += 2 * p)
            {
                uint64_t bit = (m - 1) / 2; // convert odd number m to its bit index
                // bit >> 6: bit/64
                //bit & 63: remainder for bit/64 
                // 1ULL << 6: number where only the 6th bit is a 1
                // ~(): invert. only the 6th bit is 0
                // &=: perform an AND operation: set only the 6th bit to 0; leave the others unchanged
                bitmap[bit >> 6] &= ~(1ULL << (bit & 63));
            }
        }
        if(s == 0) bitmap[0] &= ~1ULL; // 1 is not a prime number, but since it cannot be eliminated using the sieve, eliminate only the first interval directly.
        
        // count primes in this segment
        uint64_t c = 0, w;
        for(w = w0; w < w1; w++) c += __builtin_popcountll(bitmap[w]); // __builtin_popcountll is a gcc built-in function that returns the number of bits set to 1 in a 64-bit word.
        segcnt[s + 1] = c; 
    }
    double t_count = now();

    // Prefix sum of segment counts -> offsets into the prime list
    // finding the start position for writing using a cumulative sum
    for(s = 1; s <= (int64_t)nseg; s++) segcnt[s] += segcnt[s - 1];
    
    nprimes = segcnt[nseg] + 1; // segcnt[nseg] is the sum of all intervals—that is, the total number of odd prime numbers. Since the bitmap does not include 2, add 1.

    primelist = malloc(nprimes * sizeof(uint32_t));
    if(!primelist) { perror("malloc primelist"); exit(-1); }
    primelist[0] = 2;

    #pragma omp parallel for schedule(dynamic, 4)
    for(s = 0; s < (int64_t)nseg; s++)
    {
        // interval s covers the range from word w0 to, but not including, word w1
        uint64_t w0 = s * SEG_WORDS;
        uint64_t w1 = (w0 + SEG_WORDS < NWORDS) ? w0 + SEG_WORDS : NWORDS; // In the final segment, since w0 + 4096 exceeds NWORDS, so stop at NWORDS.
        uint64_t pos = 1 + segcnt[s], w; // pos = 1 + segcnt[s]: start position for writing in segment s
        for(w = w0; w < w1; w++)
        {
            uint64_t bits = bitmap[w];
            while(bits)
            {
                int t = __builtin_ctzll(bits); // position of the 1 at the very bottom
                primelist[pos++] = (uint32_t)(2 * (w * 64 + t) + 1); // Bits -> Numbers
                bits &= bits - 1; // remove the 1 at the very bottom
            }
        }
    }
    double t_end = now();

    // two largest prime numbers and semiprimes
    uint64_t P1 = primelist[nprimes - 2], P2 = primelist[nprimes - 1];
    uint64_t SP = P1 * P2;

    printf("threads=%d\n", nthreads);
    printf("start=%.9f end=%.9f total=%.6f s (base %.6f, sieve+count %.6f, list %.6f)\n", t_start, t_end, t_end - t_start, t_sieve - t_start, t_count - t_sieve, t_end - t_count);
    printf("Number of primes in [0, %llu] = %llu\n", MAX, (unsigned long long)nprimes);
    printf("Largest two primes: %llu, %llu\n", (unsigned long long)P2, (unsigned long long)P1);
    printf("SP = %llu x %llu = %llu\n", (unsigned long long)P1, (unsigned long long)P2, (unsigned long long)SP);

    // Factor test semi-primes
    // test cases 
    uint64_t tests[] = { 35ULL, 376223ULL, 4006336753ULL, 406615978649ULL, 4154092115820191ULL, 418155269059864129ULL, 999999224000135903ULL, SP }; // ULL means treat as a 64-bit unsigned integer
    int ntests = sizeof(tests) / sizeof(tests[0]), t;
    for(t = 0; t < ntests; t++)
    {
        uint64_t f1 = 0, f2 = 0;
        double t0 = now();
        int ok = factor_sp(tests[t], &f1, &f2); // pass the variable locations using &f1 and &f2, and write P1 and P2 to those locations using factor_sp
        double t1 = now();

        // if successful (ok = 1), return “SP = P1 × P2 (time)”; if unsuccessful, return “It is not the product of two prime numbers less than 1 billion.”
        if(ok) printf("SP=%llu = %llu x %llu (%.6f s)\n", (unsigned long long)tests[t], (unsigned long long)f1, (unsigned long long)f2, t1 - t0);
        else   printf("SP=%llu : not a semi-prime of primes < %llu\n", (unsigned long long)tests[t], MAX);
    }

    // Random test: pick two random primes from the list, multiply, factor
    srand(76); // 76 = random number seed to reproduce and verify the results
    for(t = 0; t < 3; t++)
    {
        // % nprimes, set the index to a value between 0 and the number of prime numbers minus 1
        uint64_t r1 = primelist[((uint64_t)rand() * 1000003ULL) % nprimes];
        uint64_t r2 = primelist[((uint64_t)rand() * 999983ULL) % nprimes];
        uint64_t f1 = 0, f2 = 0;
        int ok = factor_sp(r1 * r2, &f1, &f2); // multiply two prime numbers to form an SP, then factor it.
        printf("random: %llu x %llu = %llu -> %s (%llu x %llu)\n",
               (unsigned long long)r1, (unsigned long long)r2, (unsigned long long)(r1 * r2),
               (ok && ((f1 == r1 && f2 == r2) || (f1 == r2 && f2 == r1))) ? "PASS" : "FAIL",
               (unsigned long long)f1, (unsigned long long)f2);

        // factor_sp always puts the smaller value into f1, but sometimes the randomly selected r1 is larger
        // Therefore, even if the order is reversed, as long as they match, it’s a PASS.
    }

    // Free all memory allocated with malloc or calloc
    free(primelist); free(bitmap); free(segcnt); free(base); free(small);
    return 0;
}
