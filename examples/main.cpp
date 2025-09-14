#include <iostream>
#include "ntmath.h"
#include "utils.h"

using namespace std;
using namespace modular;

int main(/*int argc, char const *argv[]*/) {
    
    uint64_t q = 17;
    size_t N = 16;
    size_t logq = MSB(q);

    uint64_t prec = barrett_precompute(q, logq);


    uint64_t ax[N];
    uint64_t bx[N];
    uint64_t cx_add[N];
    uint64_t cx_sub[N];
    uint64_t cx_mul[N];
    uint64_t cx_mul_gold[N];

    random_uniform_vector_mod_q(ax, N, q);
    random_uniform_vector_mod_q(bx, N, q);

    cout << "q = " << q << endl;

    cout << "ax = ";
    print_array(ax, N);
    cout << endl;

    cout << "bx = ";
    print_array(bx, N);
    cout << endl;

    for (size_t i = 0; i < N; i++) {
        cx_add[i] = mod_add(ax[i], bx[i], q);
        cx_sub[i] = mod_sub(ax[i], bx[i], q);
        cx_mul[i] = mod_mul_barrett(ax[i], bx[i], q, prec, logq);

        cx_mul_gold[i] = static_cast<uint64_t>(static_cast<uint128_t>(ax[i]) * static_cast<uint128_t>(bx[i]) % static_cast<uint128_t>(q));
    }

    cout << "a + b mod q = ";
    print_array(cx_add, N);
    cout << endl;

    cout << "a - b mod q = ";
    print_array(cx_sub, N);
    cout << endl;

    cout << "a * b mod q = ";
    print_array(cx_mul, N);
    cout << endl;

    cout << "(a * b) % q = ";
    print_array(cx_mul, N);
    cout << endl;

    return 0;
}
