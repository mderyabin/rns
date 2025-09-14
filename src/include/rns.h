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

    NTL::ZZ  *Mi;     // M / mi
    uint64_t *Mi_inv; // (M / mi)^-1 mod mi

    NTL::ZZ *Bi;      // bi = Mi * mi

    void compute_M();
    void compute_Bi();
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
    NTL::ZZ from_rns(const uint64_t *rns) const;
};

}

#endif /* __RNS_H__ */