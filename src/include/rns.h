#ifndef __RNS_H__
#define __RNS_H__

#include <algorithm>
#include <NTL/ZZ.h>

#include "ntmath.h"

namespace modular {

class RNS {
private:
    size_t N;
    uint64_t *moduli;

    NTL::ZZ M;        // product of all moduli
    NTL::ZZ KMask;    // 2^K - 1, where K = bits in M

    NTL::ZZ  *Mi;     // M / mi
    uint64_t *Mi_inv; // (M / mi)^-1 mod mi
    uint64_t *Mi_inv_prec; // precomputed value for Shoup's method

    NTL::ZZ *Bi;      // bi = Mi * mi
    NTL::ZZ *Di;      // ceil(2^K * Bi / M) = ceil(2^K * Mi_inv / mi)

    size_t K;         // bits for approximate method

    void compute_M();
    void compute_Bi();
    void compute_approx();
public:

    /* main constructors and operators */
    RNS(size_t N, const uint64_t *moduli);
    ~RNS(); 

    RNS(const RNS& other);                  // Copy constructor
    RNS(RNS&& other) noexcept;              // Move constructor
    RNS& operator=(const RNS& other);       // Copy assignment operator
    RNS& operator=(RNS&& other) noexcept;   // Move assignment operator

    /* information getter */
    size_t size() const { return N; }
    const uint64_t* get_moduli() const { return moduli; }
    uint64_t get_modulus(size_t i) const { return moduli[i]; }
    NTL::ZZ get_M() const { return M; }
    size_t get_range_bits() const { return NTL::NumBits(M); }
    NTL::ZZ get_basis(size_t i) const { return Bi[i]; }

    /* conversion */
    void to_rns(uint64_t *rns, NTL::ZZ x) const;
    NTL::ZZ from_rns_classic(const uint64_t *rns) const;
    NTL::ZZ from_rns_small(const uint64_t *rns) const;
    NTL::ZZ from_rns_approx(const uint64_t *rns) const;

};

}

#endif /* __RNS_H__ */