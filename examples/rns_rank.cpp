#include <iostream>
#include "rns.h"
#include "utils.h"

using namespace std;
using namespace modular;

bool approx_double_rns_test(size_t N, size_t logq);
bool approx_int_rns_test(size_t N, size_t logq);
const size_t TEST_RNS_ITERS = 20;
const size_t TEST_NUMBER_ITERS = 100;


int main(/*int argc, char const *argv[]*/) {
    for (size_t N = 10; N < 11; N+=10) {
    
        size_t logq = 55;

        std::cout << "rank test N = " << N << ", logq = " << logq << std::endl;
        bool test_result = approx_double_rns_test(N, logq);
        if (test_result) {
            std::cout << "All tests passed!" << std::endl;
        } else {
            std::cout << "Some tests failed!" << std::endl;
        }
        std::cout << "----------------------------------------" << std::endl;
    }

    for (size_t N = 10; N < 11; N+=10) {
    
        size_t logq = 45;

        std::cout << "rank test N = " << N << ", logq = " << logq << std::endl;
        bool test_result = approx_int_rns_test(N, logq);
        if (test_result) {
            std::cout << "All tests passed!" << std::endl;
        } else {
            std::cout << "Some tests failed!" << std::endl;
        }
        std::cout << "----------------------------------------" << std::endl;
    }
}

bool approx_double_rns_test(size_t N, size_t logq) {
    uint64_t moduli[N];

    moduli[0] = FindFirstPrimeUp(logq, 2);
    for (size_t i = 1; i < N; i++) {
        moduli[i] = FindNextPrime(moduli[i-1], N);
    }

    RNS rns(N, moduli);

    bool test_passed = true;

    for (size_t i = 0; i < TEST_NUMBER_ITERS; i++) {
        // Test conversion to RNS and back
        NTL::ZZ x; // Example integer
        NTL::RandomBits(x, rns.get_range_bits() - 1); // Random integer in range [0, M)
        
        uint64_t rns_rep[N];
        rns.to_rns(rns_rep, x);

        uint64_t x_rank = rns.rank_small(rns_rep);
        uint64_t x_approx_rank = rns.rank_approx_double_small(rns_rep);
        if (x_rank != x_approx_rank) {
            // std::cout << "Test failed!" << std::endl;
            // std::cout << "x = " << x << std::endl;
            // std::cout << "x_rank = " << x_rank << std::endl;
            // std::cout << "x_approx_rank = " << x_approx_rank << std::endl;
            test_passed = false;
        } else {
            // std::cout << "Test " << i+1 << " passed." << std::endl;
        } 
    }

    return test_passed;
}

bool approx_int_rns_test(size_t N, size_t logq) {
    uint64_t moduli[N];

    moduli[0] = FindFirstPrimeUp(logq, 2);
    for (size_t i = 1; i < N; i++) {
        moduli[i] = FindNextPrime(moduli[i-1], N);
    }

    RNS rns(N, moduli);

    bool test_passed = true;

    for (size_t i = 0; i < TEST_NUMBER_ITERS; i++) {
        // Test conversion to RNS and back
        NTL::ZZ x; // Example integer
        NTL::RandomBits(x, rns.get_range_bits() - 1); // Random integer in range [0, M)
        
        uint64_t rns_rep[N];
        rns.to_rns(rns_rep, x);

        uint64_t x_rank = rns.rank_small(rns_rep);
        uint64_t x_approx_rank = rns.rank_approx_int_small(rns_rep);
        if (x_rank != x_approx_rank) {
            std::cout << "Test failed!" << std::endl;
            std::cout << "x = " << x << std::endl;
            std::cout << "x_rank = " << x_rank << std::endl;
            std::cout << "x_approx_rank = " << x_approx_rank << std::endl;
            test_passed = false;
        } else {
            std::cout << "Test " << i+1 << " passed." << std::endl;
        } 
    }

    return test_passed;
}
