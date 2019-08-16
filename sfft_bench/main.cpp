#include <iostream>
#include <utils.h>
#include "fftw3.h"
#include "sfft.h"

bool CheckAnswer(int n, const sfft_output& res, const fftw_complex* out) {
    const double EPS = 1e-3;
    for (int i = 0; i < n; ++i) {
        if (res.find(i) == res.end()) {
            if (hypot(out[i][0], out[i][1]) > EPS) {
                return false;
            }
        } else {
            if (hypot(out[i][0] - creal(res.find(i)->second), out[i][1] - cimag(res.find(i)->second)) > EPS) {
                return false;
            }
        }
    }
    return true;
}

int main() {
    int n = (1 << 13);
    int k = 50;
    srand(345);

    fftw_complex *in, *out;
    fftw_plan p;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    p = fftw_plan_dft_1d(n, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);

    sfft_plan* plan = sfft_make_plan(n, k, SFFT_VERSION_3, FFTW_ESTIMATE);
    complex_t* sin = (complex_t*) sfft_malloc(sizeof(complex_t) * n);
    {
        real_t BB = (unsigned)(1 * ((double)k));
        auto B_g1 = floor_to_pow2(BB);

        double lobefrac_g1 = 0.5 / (BB);
        int b_g1 = int (1.00 * ((double)n / B_g1));
        int w_g1;
        auto filtert1 =
            make_dolphchebyshev_t(lobefrac_g1, 1e-8, w_g1);
//        auto filterf1 =
//            make_multiple_t(filtert1, w_g1, n, b_g1).freq;
        int pos = 0;
//        while (cabs2(filterf1[pos]) > 1e-8) {
//            ++pos;
//        }
//        printf("%i\n", pos);
        printf("%i\n", w_g1);
        for (int i = 0; i < w_g1; ++i) {
            printf("%f %f\n", creal(filtert1[i]), cimag(filtert1[i]));
        }
    }
//    auto data = (sfft_v3_data*) plan->data;
//    printf("%i\n", data->w_g1);
//    for (int i = 0; i < data->w_g1; ++i) {
//        printf("%f %f\n", creal(data->filtert1[i]), cimag(data->filtert1[i]));
//    }

//    int tries = 100;
//    int ok = 0;
//    for (int iter = 0; iter < tries; ++iter) {
//
//        for (int i = 0; i < n; ++i) {
//            out[i][0] = 0;
//            out[i][1] = 0;
//        }
//
//        for (int i = 0; i < k; ++i) {
//            int pos = rand() % n;
//            out[pos][0] = 1;
////            printf("%d ", pos);
//        }
////        printf("\n");
//
//        fftw_execute(p);
//
//
//
//        for (int i = 0; i < n; ++i) {
//            sin[i] = in[i][0] + in[i][1] * I;
//        }
//
//
//        sfft_output output;
//
//        sfft_exec(plan, sin, &output);
////        printf("%ld\n", output.size());
////        for (auto w : output) {
////            printf("%d: %f %f\n", w.first, creal(w.second), cimag(w.second));
////        }
//        ok += CheckAnswer(n, output, out);
//    }
//    printf("%d/%d\n", ok, tries);
    sfft_free(sin);
//    sfft_free_plan(plan);

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
    return 0;
}