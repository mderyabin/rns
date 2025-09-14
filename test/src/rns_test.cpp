#include "gtest/gtest.h"

#include "utils.h"
#include "rns.h"

using namespace modular;

const size_t TEST_RNS_ITERS = 20;
const size_t TEST_NUMBER_ITERS = 100;

void basic_test(size_t N, size_t logq) {
        uint64_t moduli[N];

        moduli[0] = FindFirstPrimeUp(logq, 2);
        for (size_t i = 1; i < N; i++) {
            moduli[i] = FindNextPrime(moduli[i-1], N);
        }

        RNS rns(N, moduli);

        // Test size and moduli retrieval
        EXPECT_EQ(rns.size(), N);
        for (size_t i = 0; i < N; i++) {
            EXPECT_EQ(rns.get_modulus(i), moduli[i]);
            EXPECT_EQ(rns.get_basis(i) % moduli[i], 1);
        }

        for (size_t i = 0; i < TEST_NUMBER_ITERS; i++) {
            // Test conversion to RNS and back
            NTL::ZZ x; // Example integer
            NTL::RandomBits(x, rns.get_range_bits() - 1); // Random integer in range [0, M)
            
            uint64_t rns_rep[N];
            rns.to_rns(rns_rep, x);

            // Expected RNS representation: [2, 3, 2] since 23 mod 3 = 2, 23 mod 5 = 3, 23 mod 7 = 2
            for (size_t k = 0; k < N; k++) {
                EXPECT_EQ(rns_rep[k], x % moduli[k]);
            }
            NTL::ZZ x_reconstructed = rns.from_rns(rns_rep);
            EXPECT_EQ(x_reconstructed, x % rns.get_M()); // Ensure reconstruction is correct modulo M
        }
}

TEST(RNS, BasicTest) {
    for (size_t j = 0; j < TEST_RNS_ITERS; j++) {
    
        size_t N = 4 + (rand() % 21);
        size_t logq = 20 + (rand() % 35);

        basic_test(N, logq);}
}



