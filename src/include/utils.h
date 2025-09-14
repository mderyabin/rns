#ifndef __UTILS_H__
#define __UTILS_H__

#include <cstdlib>
#include <iostream>
#include <random>

namespace modular {

uint64_t generate_number(uint64_t q) {
    std::random_device rd;
    std::default_random_engine en(rd());
    std::uniform_int_distribution<uint64_t> rnd(0, q-1);

    return rnd(en);
}

void random_uniform_vector_mod_q(uint64_t *ax, size_t N, uint64_t q) {
    std::random_device rd;
    std::default_random_engine en(rd());
    std::uniform_int_distribution<uint64_t> rnd(0, q-1);

    for (size_t i = 0; i < N; i++) { 
        ax[i] = rnd(en);
    }
}

template<typename T>
void print_array(const T *ax, size_t n) {
    std::cout << "[" << ax[0];
    for (size_t i = 1; i < n; i++) {
        std::cout << ", " << ax[i];
    }
    std::cout << "]";
}

} // namespace modular 

#endif /* __UTILS_H__ */