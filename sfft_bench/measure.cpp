#include <random>
#include <cstring>
#include "iostream"
#include "chrono"
#include "fftw3.h"
#include "sfft.h"

template <class Generator>
void GenRandomSupport(int n, int d, int64_t sparsity, Generator& gen, complex_t* out) {
    int size = 1;
    for (int i= 0; i < d; ++i) {
        size *= n;
    }
    std::uniform_int_distribution<int64_t> dist(0, size - 1);
    for (int i = 0; i < sparsity; ++i) {
        out[dist(gen)] = 1;
    }
}

//std::vector<complex_t> GenDiracComb(int n, int d, int64_t sparsity) {
//    int size = 1;
//    for (int i= 0; i < d; ++i) {
//        size *= n;
//    }
//    std::uniform_int_distribution<int64_t> dist(0, size - 1);
//    std::vector<complex_t> out(size);
//    for (int i = 0; i < size; i += size / sparsity) {
//        out[i] = 1;
//    }
//    return out;
//}
//
//template <class Generator>
//std::vector<complex_t> GenCombined(const SignalInfo& info, int64_t sparsity, Generator& gen) {
//    std::vector<complex_t> out = GenRandomSupport(info, sparsity / 2, gen);
//    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / (sparsity / 2)) {
//        out[i] += 1;
//    }
//    return out;
//}

void PrintDur(const std::chrono::nanoseconds& dur) {
    std::cout << std::chrono::duration<double, std::milli>(dur).count();
}

void PrintArr(const std::string& name, const std::vector<std::chrono::nanoseconds>& durs) {
    std::cout << name << " = [\n";
    for (auto p: durs) {
        std::cout << p.count() << ",\n";
    }
    std::cout << "]\n";
}

template <class Func>
auto RunBenchmark(const std::string& name, int iters, complex_t** input, Func function, std::vector<std::chrono::nanoseconds>& res) {
    using clock = std::chrono::system_clock;
    std::cout << name << ": ";
    auto start = clock::now();
    for (int i = 0; i < iters; ++i) {
        (void) function(input[i], i);
    }
    auto dur = clock::now() - start;
    dur /= iters;
    PrintDur(dur);
    res.push_back(dur);
    std::cout << ", ";
}

int main() {
    std::mt19937_64 gen(123423);
    std::vector<int64_t> npow;
    std::vector<std::chrono::nanoseconds> dur1;

    for (int64_t pw = 6; pw <= 7; ++pw) {
        int64_t n = 1 << pw;
        int d = 3;
        int N = 1;
        for (int i = 0; i < d; ++i) {
            N *= n;
        }
        const int64_t k = 27;
        const int iters = 5;
        sfft_multidim_params params{false, 1, 0.3};

        complex_t **in, *out, *in1;
        fftw_plan p;
        out = (complex_t*) fftw_malloc(sizeof(complex_t) * N);
        in1 = (complex_t*) fftw_malloc(sizeof(complex_t) * N);
        in = (complex_t**) malloc(sizeof(complex_t*) * iters);
        std::vector<int> ranks(d, n);
        p = fftw_plan_dft(d, ranks.data(), reinterpret_cast<fftw_complex*>(out), reinterpret_cast<fftw_complex*>(in1), FFTW_BACKWARD, FFTW_ESTIMATE);
        for (int j = 0; j < iters; ++j) {
            GenRandomSupport(n, d, k, gen, out);
            fftw_execute(p);
            in[j] = (complex_t*) fftw_malloc(sizeof(complex_t) * N);
            memcpy(in[j], in1, sizeof(complex_t) * N);
        }

        sfft_plan_multidim* plan = sfft_make_plan_multidim(n, d, k, 1, params);

        sfft_output output;

        npow.push_back(pw);
        std::cout << "p = " << pw << ", ";
        RunBenchmark("sfft", iters, in, [&](complex_t* x, int){return sfft_exec_multidim(plan, x, &output), output.clear(); }, dur1);
        std::cout << std::endl;

//        sfft_free(input);
        sfft_free_plan_multidim(plan);

        fftw_destroy_plan(p);
        for (int j = 0; j < iters; ++j) {
            fftw_free(in[j]);
        }
        free(in);
        fftw_free(in1);
        fftw_free(out);
    }
    std::cout << "p = [\n";
    for (auto pw: npow) {
        std::cout << pw << ",\n";
    }
    std::cout << "]\n";
    PrintArr("sfft", dur1);

    return 0;
}

