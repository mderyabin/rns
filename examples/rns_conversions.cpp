#include <iostream>
#include "rns.h"
#include "utils.h"

using namespace std;
using namespace modular;

int main(/*int argc, char const *argv[]*/) {
    
    const size_t N = 8;
    const size_t logq = 50;

    uint64_t moduli[N];

    moduli[0] = FindFirstPrimeDown(logq, N);
    for (size_t i = 1; i < N; i++) {
        moduli[i] = FindPrevPrime(moduli[i-1], N);
    }

    cout << "moduli = ";
    print_array(moduli, N);
    cout << endl;

    RNS rns(N, moduli);

    NTL::ZZ X;
    NTL::RandomBits(X, (N-1)*logq); 

    cout << "x = " << X << endl;

    uint64_t rns_x[N];
    rns.to_rns(rns_x, X);

    cout << "x in RNS = ";
    print_array(rns_x, N);
    cout << endl;

    NTL::ZZ x_recovered = rns.from_rns_classic(rns_x);
    cout << "x recovered from RNS = " << x_recovered << endl;

    return 0;
}