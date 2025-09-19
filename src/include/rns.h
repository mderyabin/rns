#ifndef __RNS_H__
#define __RNS_H__

#include <algorithm>
#include <NTL/ZZ.h>

#include "ntmath.h"

namespace modular {

class RNS {
public:
    size_t N;
    uint64_t *moduli;

    NTL::ZZ M;        // product of all moduli
    NTL::ZZ KMask;    // 2^K - 1, where K = bits in M

    NTL::ZZ  *Mi;     // M / mi
    uint64_t *Mi_inv; // (M / mi)^-1 mod mi
    uint64_t *Mi_inv_prec; // precomputed value for Shoup's method

    NTL::ZZ *Bi;      // bi = Mi * mi
    NTL::ZZ *Di;      // ceil(2^K * Bi / M) = ceil(2^K * Mi_inv / mi)

    double *m_inv_approx;           // precomputed values for approx double method to compute rank
    double *m_approx;               // precomputed values for approx double method to compute rank
    uint64_t *m_inv_approx_int;     // precomputed values for approx int method to compute rank

    size_t K;         // bits for approximate method
    size_t K_rank;    // bits for approximate rank method

    void compute_M();
    void compute_Bi();
    void compute_approx();
    void compute_rank();
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

    /* rank */
    uint64_t rank_small(const uint64_t *rns) const;
    uint64_t rank_approx_double_small(const uint64_t *rns) const;
    uint64_t rank_approx_int_small(const uint64_t *rns) const;

};

}

#endif /* __RNS_H__ */