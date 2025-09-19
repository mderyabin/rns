#include "benchmark/benchmark.h"

#include "rns.h"
#include "utils.h"

using namespace modular;

void BC_RNS_Small_Rank(benchmark::State& state, size_t N, size_t logq) {

    uint64_t moduli[N];

    moduli[0] = FindFirstPrimeDown(logq, 1<<18);
    for (size_t i = 1; i < N; i++) {
        moduli[i] = FindPrevPrime(moduli[i-1], 1<<18);
    }

    RNS rns(N, moduli);

    NTL::ZZ X;
    NTL::RandomBits(X, rns.get_range_bits()); 

    uint64_t rns_x[N];
    rns.to_rns(rns_x, X);

    uint64_t x_recovered;

    while (state.KeepRunning()) {
        benchmark::DoNotOptimize(rns);

        for (size_t i = 0; i < N; ++i) {
            benchmark::DoNotOptimize(rns_x[i]);
        }

        x_recovered = rns.rank_small(rns_x);

        benchmark::DoNotOptimize(x_recovered);
    }
}

void BC_RNS_ApproxDouble_Small_Rank(benchmark::State& state, size_t N, size_t logq) {

    uint64_t moduli[N];

    moduli[0] = FindFirstPrimeDown(logq, 1<<18);
    for (size_t i = 1; i < N; i++) {
        moduli[i] = FindPrevPrime(moduli[i-1], 1<<18);
    }

    RNS rns(N, moduli);

    NTL::ZZ X;
    NTL::RandomBits(X, rns.get_range_bits()); 

    uint64_t rns_x[N];
    rns.to_rns(rns_x, X);

    uint64_t x_recovered;

    while (state.KeepRunning()) {
        benchmark::DoNotOptimize(rns);

        benchmark::DoNotOptimize(rns.Mi_inv);
        benchmark::DoNotOptimize(rns.Mi_inv_prec);
        benchmark::DoNotOptimize(rns.moduli);

        for (size_t i = 0; i < N; ++i) {
            benchmark::DoNotOptimize(rns_x[i]);
        }
        
        x_recovered = rns.rank_approx_double_small(rns_x);

        benchmark::DoNotOptimize(x_recovered);
    }
}

void BC_RNS_ApproxInt_Small_Rank(benchmark::State& state, size_t N, size_t logq) {

    uint64_t moduli[N];

    moduli[0] = FindFirstPrimeDown(logq, 1<<18);
    for (size_t i = 1; i < N; i++) {
        moduli[i] = FindPrevPrime(moduli[i-1], 1<<18);
    }

    RNS rns(N, moduli);

    NTL::ZZ X;
    NTL::RandomBits(X, rns.get_range_bits()); 

    uint64_t rns_x[N];
    rns.to_rns(rns_x, X);

    uint64_t x_recovered;

    while (state.KeepRunning()) {
    
        
        benchmark::DoNotOptimize(rns);
        benchmark::DoNotOptimize(rns.Mi_inv);
        benchmark::DoNotOptimize(rns.Mi_inv_prec);
        benchmark::DoNotOptimize(rns.moduli);


        for (size_t i = 0; i < N; ++i) {
            benchmark::DoNotOptimize(rns_x[i]);
        }
        
        x_recovered = rns.rank_approx_int_small(rns_x);

        benchmark::DoNotOptimize(x_recovered);
    }
}

// BENCHMARK_CAPTURE(BC_RNS_Small_Rank, Small_10_50,   10, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
// BENCHMARK_CAPTURE(BC_RNS_Small_Rank, Small_20_50,   20, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
// BENCHMARK_CAPTURE(BC_RNS_Small_Rank, Small_30_50,   30, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
// BENCHMARK_CAPTURE(BC_RNS_Small_Rank, Small_40_50,   40, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
// BENCHMARK_CAPTURE(BC_RNS_Small_Rank, Small_50_50,   50, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
// BENCHMARK_CAPTURE(BC_RNS_Small_Rank, Small_60_50,   60, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
// BENCHMARK_CAPTURE(BC_RNS_Small_Rank, Small_70_50,   70, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
// BENCHMARK_CAPTURE(BC_RNS_Small_Rank, Small_80_50,   80, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
// BENCHMARK_CAPTURE(BC_RNS_Small_Rank, Small_90_50,   90, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
// BENCHMARK_CAPTURE(BC_RNS_Small_Rank, Small_100_50, 100, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);

BENCHMARK_CAPTURE(BC_RNS_ApproxDouble_Small_Rank, ApproxDouble_10_50,   10, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_ApproxDouble_Small_Rank, ApproxDouble_20_50,   20, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_ApproxDouble_Small_Rank, ApproxDouble_30_50,   30, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_ApproxDouble_Small_Rank, ApproxDouble_40_50,   40, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_ApproxDouble_Small_Rank, ApproxDouble_50_50,   50, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_ApproxDouble_Small_Rank, ApproxDouble_60_50,   60, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_ApproxDouble_Small_Rank, ApproxDouble_70_50,   70, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_ApproxDouble_Small_Rank, ApproxDouble_80_50,   80, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_ApproxDouble_Small_Rank, ApproxDouble_90_50,   90, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_ApproxDouble_Small_Rank, ApproxDouble_100_50, 100, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);

BENCHMARK_CAPTURE(BC_RNS_ApproxInt_Small_Rank, ApproxInt_10_50,   10, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_ApproxInt_Small_Rank, ApproxInt_20_50,   20, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_ApproxInt_Small_Rank, ApproxInt_30_50,   30, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_ApproxInt_Small_Rank, ApproxInt_40_50,   40, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_ApproxInt_Small_Rank, ApproxInt_50_50,   50, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_ApproxInt_Small_Rank, ApproxInt_60_50,   60, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_ApproxInt_Small_Rank, ApproxInt_70_50,   70, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_ApproxInt_Small_Rank, ApproxInt_80_50,   80, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_ApproxInt_Small_Rank, ApproxInt_90_50,   90, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_ApproxInt_Small_Rank, ApproxInt_100_50, 100, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);


BENCHMARK_MAIN();