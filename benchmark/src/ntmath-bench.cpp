#include "benchmark/benchmark.h"

#include "ntmath.h"
#include "utils.h"

using namespace modular;

void BC_ModAddNaive(benchmark::State& state) {
    uint64_t q = generate_number(1ULL<<63);

    uint64_t a = generate_number(q);

    uint64_t b = generate_number(q);

    uint64_t c;


    while (state.KeepRunning()) {
        benchmark::DoNotOptimize(a);
        benchmark::DoNotOptimize(b);

        c = (a+b) % q;

        benchmark::DoNotOptimize(c);
    }
}

void BC_ModAdd(benchmark::State& state) {
    uint64_t q = generate_number(1ULL<<62);

    uint64_t a = generate_number(q);

    uint64_t b = generate_number(q);

    uint64_t c;


    while (state.KeepRunning()) {
        benchmark::DoNotOptimize(a);
        benchmark::DoNotOptimize(b);
        benchmark::DoNotOptimize(q);
        
        c = mod_add(a, b, q);

        benchmark::DoNotOptimize(c);
    }
}

void BC_ModSubNaive(benchmark::State& state) {
    uint64_t q = generate_number(1ULL<<62);

    uint64_t a = generate_number(q);

    uint64_t b = generate_number(q);

    uint64_t c;


    while (state.KeepRunning()) {
        benchmark::DoNotOptimize(a);
        benchmark::DoNotOptimize(b);
        benchmark::DoNotOptimize(q);

        c = (a-b) % q;

        benchmark::DoNotOptimize(c);
    }
}

void BC_ModSub(benchmark::State& state) {
    uint64_t q = generate_number(1ULL<<62);

    uint64_t a = generate_number(q);

    uint64_t b = generate_number(q);

    uint64_t c;


    while (state.KeepRunning()) {
        benchmark::DoNotOptimize(a);
        benchmark::DoNotOptimize(b);
        benchmark::DoNotOptimize(q);

        c = mod_sub(a, b, q);

        benchmark::DoNotOptimize(c);
    }
}

void BC_ModMulNaive(benchmark::State& state) {
    uint64_t q = generate_number(1ULL<<62);

    uint64_t a = generate_number(q);

    uint64_t b = generate_number(q);

    uint64_t c;


    while (state.KeepRunning()) {
        benchmark::DoNotOptimize(a);
        benchmark::DoNotOptimize(b);
        benchmark::DoNotOptimize(q);

        c = static_cast<uint64_t>(static_cast<uint128_t>(a) * static_cast<uint128_t>(b) % static_cast<uint128_t>(q));

        benchmark::DoNotOptimize(c);
    }
}

void BC_ModMulBarrett(benchmark::State& state) {
    uint64_t q = generate_number(1ULL<<62);
    size_t logq = MSB(q);
    uint64_t prec = barrett_precompute(q, logq);

    uint64_t a = generate_number(q);

    uint64_t b = generate_number(q);

    uint64_t c;


    while (state.KeepRunning()) {
        benchmark::DoNotOptimize(a);
        benchmark::DoNotOptimize(b);
        benchmark::DoNotOptimize(q);
        benchmark::DoNotOptimize(prec);
        benchmark::DoNotOptimize(logq);

        c = mod_mul_barrett(a, b, q, prec, logq);

        benchmark::DoNotOptimize(c);
    }
}

void BC_ModMulShoup(benchmark::State& state) {
    uint64_t q = generate_number(1ULL<<62);

    uint64_t a = generate_number(q);

    uint64_t b = generate_number(q);
    uint64_t b_prec = shoup_precompute(b, q);

    uint64_t c;


    while (state.KeepRunning()) {
        benchmark::DoNotOptimize(a);
        benchmark::DoNotOptimize(b);
        benchmark::DoNotOptimize(q);
        benchmark::DoNotOptimize(b_prec);

        c = mod_mul_shoup(a, b, q, b_prec);

        benchmark::DoNotOptimize(c);
    }
}

BENCHMARK(BC_ModAddNaive)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK(BC_ModAdd)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK(BC_ModSubNaive)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK(BC_ModSub)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);

BENCHMARK(BC_ModMulNaive)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK(BC_ModMulBarrett)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);
BENCHMARK(BC_ModMulShoup)->Unit(benchmark::kNanosecond)->Repetitions(10)->ReportAggregatesOnly(true);



BENCHMARK_MAIN();