/*
Copyright 2016 Andrew R. Booker

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>

typedef int32_t i32;
typedef int64_t i64;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;
typedef __uint128_t u128;

static i64 *psptable63;
static i32 *psptable63_index;
static u32 *collision_table;
static i32 *collision_index;
static u16 *factor_table;
static u32 smallprime[]={3,5,7,11,13,17,19,23,29,31,37,41,43,47,53};

const u16 *initfactor63(const char *db) {
	int fd;
	u64 length;
	void *map;

	if ((fd=open(db,O_RDONLY)) < 0)
		return (u16 *)0;
	length = 23355139ul*sizeof(i64)+64ul*sizeof(i32)+146144317ul*sizeof(u32)
		+808315ul*sizeof(i32)+1518500250ul*sizeof(u16);
	map =  mmap(NULL,length,PROT_READ,MAP_SHARED,fd,0);
	if (map == MAP_FAILED)
		return (u16 *)0;
	madvise(map,length,MADV_WILLNEED);

	psptable63 = (i64 *)map;
	psptable63_index = (i32 *)&psptable63[23355139];
	collision_table = (u32 *)&psptable63_index[64];
	collision_index = (i32 *)&collision_table[146144317];
	factor_table = (u16 *)&collision_index[808315];

	return factor_table;
}

static inline uint64_t _uint128_mulredc(uint64_t x, uint64_t y, uint64_t n, uint64_t npi)
{
  union { u128 d; u64 l[2]; } t;
  t.d = (u128) x*y;
  t.d += (t.l[0] * npi) * (u128) n;
  t.l[1] -= n;
  return t.l[1] + (n & ((i64)t.l[1] >> 63));
}
#define mulredc(x,y)  _uint128_mulredc(x,y,n,nbar)

// returns true if n passes a strong Fermat test to base 2
// assumes n is odd and 1 < n < 2^63
static inline int isprobableprime(u64 n) {
	u64 one,neg1,x,nbar;
	i32 k;
	i64 q;

	// compute nbar = -n^{-1} mod 2^64 by Newton's method
	k = ((n+1)>>2<<3)-n; // correct mod 2^4
	k *= 2+k*(i32)n; k *= 2+k*(i32)n; k *= 2+k*(i32)n;
	nbar = (u64)k*(2+(u64)k*n);

	// compute Montgomery representations of 1 and -1
	one = 1+(~0ul)%n, neg1 = n-one;

	q = n>>1; k = __builtin_ctzl(q);
	q <<= __builtin_clzl(q);
	x = one-neg1; x += n&((i64)x>>63);
	for (q<<=1;q;q<<=1) {
		x = mulredc(x,x);
		if (q < 0) { x += x-n; x += n&((i64)x>>63); }
	}
	if (x == one || x == neg1) return 1;
	while (--k >= 0) {
		x = mulredc(x,x);
		if (x == neg1) return 1;
	}
	return 0;
}

// skips trial division
// n must be odd and positive
int fastisprime63(i64 n) {
	i64 i,j,t;

	if (!isprobableprime((u64)n)) return 0;
	t = 63-__builtin_clzl(n);
	i = psptable63_index[t]; j = psptable63_index[t+1];
	while (j-i > 1) {
		t = (i+j)>>1;
		if (psptable63[t] > n) j = t; else i = t;
	}
	return (psptable63[i] != n);
}

#define M 3037000500u
int isprime63(i64 n) {
	int i;

	if (n <= 1) return 0;
	if (!(n&1)) return (n==2);
	if (n < M) return !factor_table[n>>1];
	for (i=0;i<sizeof(smallprime)/sizeof(smallprime[0]);i++)
		if (!(n % smallprime[i]))
			return (n == smallprime[i]);
	return fastisprime63(n);
}

// given that f is odd, less than M and divides n,
// return a list of the prime factors of f
// and remove them from n
static int smallfactors63(i64 *p,int *e,u32 f,u64 *n) {
	u32 y;
	int k;

	k = 0;
	if (f > 1) {
		*n /= f;
		do {
			if (!(y=factor_table[f>>1])) y = f;
			p[k] = y, e[k] = 0;
			do f /= y, e[k]++; while (f % y == 0);
			while (*n % y == 0)
				*n /= y, e[k]++;
			k++;
		} while (f > 1);
	}
	return k;
}

// assumes x is odd
static inline u64 oddgcd(u64 x,u64 y) {
	if (y) {
		y >>= __builtin_ctzl(y);
		while (x != y)
			if (x < y)
				y -= x, y >>= __builtin_ctzl(y);
			else
				x -= y, x >>= __builtin_ctzl(x);
	}
	return x;
}

int factor63(i64 *p,int *e,i64 n0) {
// iterations can be any number > 13193 such that
//   iterations % maxstride == 0 && collision_index[iterations-1] != 0
// maxstride must be a power of 2
#define iterations 300288
#define maxstride 256
	u64 m,n,nbar,x,y,y0,f,one;
	u32 *ptr,s;
	int i,j,k,mask;

	k = 0;
	if (n0 < 0) {
		p[k] = -1, e[k] = 1;
		k++;
		n = -n0;
	} else
		n = n0;

	if (!(n & 1)) {
		p[k] = 2, e[k] = __builtin_ctzl(n);
		n >>= e[k];
		k++;
	}
	f = oddgcd(n,16294579238595022365ul);
	for (i=0;f>1;i++)
		if (f % smallprime[i] == 0) {
			f /= smallprime[i];
			p[k] = smallprime[i], e[k] = 0;
			do n /= smallprime[i], e[k]++; while (n % smallprime[i] == 0);
			k++;
		}
	if (n < M) {
		k += smallfactors63(p+k,e+k,n,&n);
		return k;
	}
	if (fastisprime63(n)) {
		p[k] = n, e[k] = 1;
		k++;
		return k;
	}

	// compute nbar = -n^{-1} mod 2^64 by Newton's method
	i = ((n+1)>>2<<3)-n; // correct mod 2^4
	i *= 2+i*(int)n; i *= 2+i*(int)n; i *= 2+i*(int)n;
	nbar = (u64)i*(2+(u64)i*n);

	// compute Montgomery representation of 1
	one = 1+(~0ul)%n;

	// Pollard rho
	m = n, y = f = one;
	for (i=1;i<iterations;i<<=1) {
		mask = (i < maxstride) ? i-1 : maxstride-1;
		x = y0 = y; j = 0;
		do {
			y = mulredc(y,y+one);
			f = mulredc(f,labs(y-x));
			j++;
			if (j & mask) continue;
			if ((f=oddgcd(m,f))==1) {
				y0 = y;
				continue;
			}
			if (f >= M) {
				// backtrack to find exact cycle time
				y = y0, j -= (mask+1);
				do {
					y = mulredc(y,y+one);
					f = oddgcd(m,labs(y-x));
					j++;
				} while (f == 1);

				ptr = collision_table+collision_index[i+j-2];
				while (f >= M)
					if (fastisprime63(f)) {
						p[k] = f, e[k] = 1;
						m /= f; k++; f = 1;
					} else if ((s=(u32)sqrtl((long double)f)),f==(u64)s*s)
						f = s;
					else {
						while (f % *ptr) ptr++;
						p[k] = *ptr++; e[k] = 0;
						do f /= p[k], m /= p[k], e[k]++; while (f % p[k] == 0);
						while (m % p[k] == 0) m /= p[k], e[k]++;
						k++;
					}
			}
			k += smallfactors63(p+k,e+k,f,&m);
			if (m < M) {
				k += smallfactors63(p+k,e+k,m,&m);
				return k;
			}
			if (fastisprime63(m)) {
				p[k] = m, e[k] = 1;
				k++;
				return k;
			}

			if (!(j&mask)) y0 = y;
			f = one;
		} while (j < i && i+j < iterations);
	}

	// the only remaining primes are those with long cycles
	// since iterations > 13193, every prime below 2^21 has been considered,
	//   so m is either a prime square or semiprime
	ptr = collision_table+collision_index[iterations-1];
	while (m % *ptr) ptr++;
	p[k] = *ptr; m /= p[k];
	if (m == p[k])
		e[k] = 2;
	else {
		e[k++] = 1;
		p[k] = m, e[k] = 1;
	}
	k++;
	return k;
}

#if 0
#include <stdio.h>

static void sort_factors(int n, i64 *f, int *e) {
  int i, j, te;
  i64 tf;
  for (i = 1; i < n; i++)
    for (j = i; j > 0 && f[j-1] > f[j]; j--)
      { tf=f[j-1]; f[j-1]=f[j]; f[j]=tf;  te=e[j-1]; e[j-1]=e[j]; e[j]=te; }
}

int main(int argc,char *argv[]) {
	i64 n,p[16];
	int e[16];
	int i,j,k;

	if (!initfactor63("factor.bin")) {
		fprintf(stderr, "Cannot read factor data\n");
		return -1;
	}

	if (argc <= 1) {  /* no args, read from stdin */
		char line[1024];
		while (fgets(line, sizeof(line), stdin)) {
			unsigned long n = strtoul(line, 0, 10);
			int nfactors = factor63(p, e, (i64) n);
			sort_factors(nfactors, p, e);
			printf("%lu:", n);
			for (i = 0; i < nfactors; i++)
				for (j = 0; j < e[i]; j++)
					printf(" %lld", p[i]);
			printf("\n");
		}
		return(0);
	}

	for (i=1;i<argc;i++)
		if (sscanf(argv[i],"%ld",(long*)&n) == 1) {
			k = factor63(p,e,n);
			printf("%ld",(long)n);
			for (j=0;j<k;j++) {
				printf(" %c %ld",j?'*':'=',(long)p[j]);
				if (e[j] > 1) printf("^%d",e[j]);
			}
			printf("\n");
		}

	return 0;
}
#endif
