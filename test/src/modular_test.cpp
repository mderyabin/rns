#include "gtest/gtest.h"

#include "utils.h"
#include "ntmath.h"

using namespace modular;

const size_t N = 10000;
const size_t logq = 58;

TEST(ModularArithmetic, ModAddTest) {
    uint64_t q = generate_number(1ULL<<logq);
    
    uint64_t ax[N];
    uint64_t bx[N];
    uint64_t cx[N];
    uint64_t cx_gold[N];

    random_uniform_vector_mod_q(ax, N, q);
    random_uniform_vector_mod_q(bx, N, q);

    /* testing add */
    for (size_t i = 0; i < N; i++) {
        cx[i] = mod_add(ax[i], bx[i], q);
        cx_gold[i] = (ax[i] + bx[i]) % q;

        EXPECT_EQ(cx[i], cx_gold[i]);
    }
}

TEST(ModularArithmetic, ModSubTest) {
    uint64_t q = generate_number(1ULL<<logq);
    
    uint64_t ax[N];
    uint64_t bx[N];
    uint64_t cx[N];
    uint64_t cx_gold[N];

    random_uniform_vector_mod_q(ax, N, q);
    random_uniform_vector_mod_q(bx, N, q);

    /* testing sub */
    for (size_t i = 0; i < N; i++) {
        cx[i] = mod_sub(ax[i], bx[i], q);
        cx_gold[i] = (q + ax[i] - bx[i]) % q;

        EXPECT_EQ(cx[i], cx_gold[i]);
    }
}

TEST(ModularArithmetic, ModBarrettTest) {
    uint64_t q = generate_number(1ULL<<logq);
    
    uint64_t ax[N];
    uint64_t bx[N];
    uint64_t cx[N];
    uint64_t cx_gold[N];

    random_uniform_vector_mod_q(ax, N, q);
    random_uniform_vector_mod_q(bx, N, q);

    /* testing barrett */
    uint64_t prec = barrett_precompute(q, logq);

    for (size_t i = 0; i < N; i++) {
        cx[i] = mod_mul_barrett(ax[i], bx[i], q, prec, logq);
        cx_gold[i] = static_cast<uint64_t>(static_cast<uint128_t>(ax[i]) * static_cast<uint128_t>(bx[i]) % static_cast<uint128_t>(q));    
    
        EXPECT_EQ(cx[i], cx_gold[i]);
    }
}

TEST(ModularArithmetic, ModShoupTest) {
    uint64_t q = generate_number(1ULL<<logq);

    uint64_t ax[N];
    uint64_t bx[N];
    uint64_t cx[N];
    uint64_t cx_gold[N];

    random_uniform_vector_mod_q(ax, N, q);
    random_uniform_vector_mod_q(bx, N, q);

    /* testing shoup */
    for (size_t i = 0; i < N; i++) {
        uint64_t b_prec = shoup_precompute(bx[i], q);

        cx[i] = mod_mul_shoup(ax[i], bx[i], q, b_prec);
        cx_gold[i] = static_cast<uint64_t>(static_cast<uint128_t>(ax[i]) * static_cast<uint128_t>(bx[i]) % static_cast<uint128_t>(q));
    
        EXPECT_EQ(cx[i], cx_gold[i]);
    }
}