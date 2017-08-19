This is the documentation for factor63, a library for primality testing
and factoring of integers less than 2^63, by Andrew R. Booker. It is free
software, distributed under the terms of the GPLv3. For more information,
see the file LICENSE.txt.

Prior to using the library, you must download and decompress the
associated data table from
http://people.maths.bris.ac.uk/~maarb/code/factor.bin.xz

The source code for the library consists of a single C file, factor63.c.
It exports four functions:

const uint16_t *initfactor63(const char *filename);
Initializes the data table pointed to by filename. The table uses about
3.8GB of memory, but the file is accessed via memory mapping, so separate
processes/threads will share a single copy.

The return value is NULL if there is an error. Otherwise it points
to a table that is useful for primality testing and factoring small
integers, even without the library. For any odd positive integer
n < 3037000500, table[(n-1)/2] is defined to be 0 if n is prime or 1,
and the smallest prime factor of n otherwise.

int isprime63(int64_t n);
Returns 1 if n is positive and prime, and 0 otherwise.

int fastisprime63(int64_t n);
Version of the above that skips the positivity check and trial division.
The input n must be odd and positive (results are undefined otherwise).

int factor63(int64_t *p,int *e,int64_t n);
Places the prime factors of n in p and their exponents in e, and returns
the number of factors. If n is negative, p[0] is set to -1 and e[0] is
set to 1.
