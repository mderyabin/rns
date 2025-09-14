#include "rns.h"

using namespace NTL;

namespace modular {

void RNS::compute_M() {
    M = ZZ(1);
    for (size_t i = 0; i < N; i++) {
        M *= moduli[i];
    }
}
void RNS::compute_Bi() {
    Mi = new ZZ[N];
    Mi_inv = new uint64_t[N];
    Bi = new ZZ[N];
    for (size_t i = 0; i < N; i++) {
        Mi[i] = M / moduli[i];
        Mi_inv[i] = InvMod(Mi[i] % moduli[i], moduli[i]);
        Bi[i] = Mi[i] * Mi_inv[i];
    }
}

RNS::RNS(size_t N, const uint64_t *moduli) : N(N) {
    this->moduli = new uint64_t[N];
    for (size_t i = 0; i < N; i++) {
        this->moduli[i] = moduli[i];
    }

    compute_M();
    compute_Bi();
}

RNS::~RNS() {
    delete[] moduli;
    delete[] Mi;
    delete[] Mi_inv;
    delete[] Bi;
}

// Copy constructor
RNS::RNS(const RNS& other) : N(other.N) {
    if (other.moduli) {
        moduli = new uint64_t[N];
        std::copy(other.moduli, other.moduli + N, moduli);
    } else {
        moduli = nullptr;
    }
}

// Move constructor
RNS::RNS(RNS&& other) noexcept : N(other.N), moduli(other.moduli) {
    other.N = 0;
    other.moduli = nullptr;
}

// Copy assignment operator
RNS& RNS::operator=(const RNS& other) {
    if (this != &other) {
        delete[] moduli;
        N = other.N;
        if (other.moduli) {
            moduli = new uint64_t[N];
            std::copy(other.moduli, other.moduli + N, moduli);
        } else {
            moduli = nullptr;
        }
    }
    return *this;
}

// Move assignment operator
RNS& RNS::operator=(RNS&& other) noexcept {
    if (this != &other) {
        delete[] moduli;
        N = other.N;
        moduli = other.moduli;
        other.N = 0;
        other.moduli = nullptr;
    }
    return *this;
}

// simple conversion functions to RNS from big integer
// don't have to be efficient
void RNS::to_rns(uint64_t *rns, NTL::ZZ x) const {
    NTL::ZZ r;
    for (size_t i = 0; i < N; i++) {
        r = x % moduli[i];
        // rns[i] = conv<unsigned long long>(r);
        long t;
        conv(t, r);
        rns[i] = static_cast<uint64_t>(t);
    }
}

// simple conversion functions from RNS to big integer
// don't have to be efficient
NTL::ZZ RNS::from_rns(const uint64_t *rns) const {
    NTL::ZZ x = ZZ(0);
    for (size_t i = 0; i < N; i++) {
        x += Bi[i] * rns[i]; 
    }
    return x % M;
}

}