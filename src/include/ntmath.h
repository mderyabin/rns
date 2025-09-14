#ifndef __NTMATH_H__
#define __NTMATH_H__

#include <cstdlib>
#include <vector>
#include <unordered_set>
#include <exception>

namespace modular {

typedef unsigned __int128 uint128_t;
size_t MSB(uint64_t t);

uint64_t barrett_precompute(uint64_t q, size_t logq);
uint64_t mod_mul_barrett(uint64_t a, uint64_t b, uint64_t q, uint64_t prec, size_t logq);

uint64_t shoup_precompute(uint64_t c, uint64_t q);
uint64_t mod_mul_shoup(uint64_t a, uint64_t c, uint64_t q, uint64_t prec);

// c = a+b mod q
inline uint64_t mod_add(uint64_t a, uint64_t b, uint64_t q) {
    // uint64_t d = a + b;
    // return (d < q) ? d : d - q;
    return (a + b < q) ? a + b : a + b - q;
}

// c = a-b mod q
inline uint64_t mod_sub(uint64_t a, uint64_t b, uint64_t q) {
    // uint64_t d = a - b;
    // return (a > b) ? d : q + d;
    return (a >= b) ? a - b : q + a - b;
}


/* semi-efficient helper functions */

uint64_t mod_exp(uint64_t a, uint64_t e, uint64_t m, uint64_t prec = 0, size_t logm = 0);

bool IsPrime(uint64_t m, size_t iters = 1000); // Miller-Rabin primarity test

uint64_t FindFirstPrimeDown(size_t logm, size_t M); // find prime in form m = (2^logm + 1) - k*M
uint64_t FindFirstPrimeUp(size_t logm, size_t M); // find prime in form m = (2^logm + 1) - k*M

uint64_t FindPrevPrime(uint64_t m, size_t M); // find previous prime in form m1 = m - k*M
uint64_t FindNextPrime(uint64_t m, size_t M); // find next prime in form m = m + k*M

uint64_t GCD(uint64_t a, uint64_t b);

uint64_t RhoPollard(uint64_t a);

std::vector<uint64_t> Factorize(uint64_t a);

uint64_t FindPrimitive(uint64_t n);
uint64_t FindGenerator(uint64_t m, size_t M);

uint64_t ModInvPrime(uint64_t a, uint64_t m);



}

#endif /* __NTMATH_H__ */