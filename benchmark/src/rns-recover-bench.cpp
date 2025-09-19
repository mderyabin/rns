#include "benchmark/benchmark.h"

#include "rns.h"
#include "utils.h"

using namespace modular;


void BC_RNS_Classic_Convert(benchmark::State& state, size_t N, size_t logq) {

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

    NTL::ZZ x_recovered;

    while (state.KeepRunning()) {
        benchmark::DoNotOptimize(rns);
        for (size_t i = 0; i < N; ++i) {
            benchmark::DoNotOptimize(rns_x[i]);
        }

        x_recovered = rns.from_rns_classic(rns_x);

        benchmark::DoNotOptimize(x_recovered);
    }
}

void BC_RNS_Small_Convert(benchmark::State& state, size_t N, size_t logq) {

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

    NTL::ZZ x_recovered;

    while (state.KeepRunning()) {
        benchmark::DoNotOptimize(rns);

        for (size_t i = 0; i < N; ++i) {
            benchmark::DoNotOptimize(rns_x[i]);
        }

        x_recovered = rns.from_rns_small(rns_x);

        benchmark::DoNotOptimize(x_recovered);
    }
}

void BC_RNS_Approx_Convert(benchmark::State& state, size_t N, size_t logq) {

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

    NTL::ZZ x_recovered;

    while (state.KeepRunning()) {
        benchmark::DoNotOptimize(rns);

        for (size_t i = 0; i < N; ++i) {
            benchmark::DoNotOptimize(rns_x[i]);
        }

        x_recovered = rns.from_rns_approx(rns_x);

        benchmark::DoNotOptimize(x_recovered);
    }
}


BENCHMARK_CAPTURE(BC_RNS_Classic_Convert, Classic_10_50, 10, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Classic_Convert, Classic_20_50, 20, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Classic_Convert, Classic_30_50, 30, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Classic_Convert, Classic_40_50, 40, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Classic_Convert, Classic_50_50, 50, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Classic_Convert, Classic_60_50, 60, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Classic_Convert, Classic_70_50, 70, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Classic_Convert, Classic_80_50, 80, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Classic_Convert, Classic_90_50, 90, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Classic_Convert, Classic_100_50, 100, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);

BENCHMARK_CAPTURE(BC_RNS_Small_Convert, Small_10_50, 10, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Small_Convert, Small_20_50, 20, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Small_Convert, Small_30_50, 30, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Small_Convert, Small_40_50, 40, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Small_Convert, Small_50_50, 50, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Small_Convert, Small_60_50, 60, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Small_Convert, Small_70_50, 70, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Small_Convert, Small_80_50, 80, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Small_Convert, Small_90_50, 90, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Small_Convert, Small_100_50, 100, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);

BENCHMARK_CAPTURE(BC_RNS_Approx_Convert, Approx_10_50, 10, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Approx_Convert, Approx_20_50, 20, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Approx_Convert, Approx_30_50, 30, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Approx_Convert, Approx_40_50, 40, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Approx_Convert, Approx_50_50, 50, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Approx_Convert, Approx_60_50, 60, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Approx_Convert, Approx_70_50, 70, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Approx_Convert, Approx_80_50, 80, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Approx_Convert, Approx_90_50, 90, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK_CAPTURE(BC_RNS_Approx_Convert, Approx_100_50, 100, 50)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);

BENCHMARK_MAIN();