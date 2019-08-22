#include <iostream>
#include <utils.h>
#include "fftw3.h"
#include "sfft.h"
#include "fftw.h"
#include "vector"

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
    const int n = (1 << 7);
    const int d = 3;
    int N = 1;
    for (int i = 0; i < d; ++i) {
        N *= n;
    }
    int k = 27;
    srand(8947);
    fftw_complex *in, *out;
    fftw_plan p;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    std::vector<int> ranks(d, n);
    p = fftw_plan_dft(d, ranks.data(), out, in, FFTW_BACKWARD, FFTW_ESTIMATE);

    sfft_plan_multidim* plan = sfft_make_plan_multidim(n, d, k, 1);
    complex_t* input = (complex_t*) sfft_malloc(sizeof(complex_t) * N);

//    printf("%d\n", plan->data.filters[0].B_g);
//    printf("%d\n", plan->data.filters[0].sizef);
//    printf("%d\n", plan->data.filters[0].sizet);
//    printf("%d %d\n", n / plan->data.filters[0].B_g, n);
//    int pos = 0;
//    while (cabs2(plan->data.filters[0].freq_at(pos) - 1) < 1e-3) {
//        ++pos;
//    }
//    printf("%d\n", pos * 2);

//    complex_t* x = (complex_t*) fftw_malloc(n * sizeof(complex_t));
//    for (int i = 0; i < n; ++i) {
//        printf("%f %f\n", creal(plan->data.filters[0].freq_at(i)), cimag(plan->data.filters[0].freq_at(i)));
//        x[i] = plan->data.filters[0].time_at(i);
//    }

//    fftw_dft(x, n, x);
//    int cnt = 0;
//    for (int i = 0; i < n; ++i) {
//        cnt += (cabs2(x[i] - plan->data.filters[0].freq_at(i)) < 1e-6);
//        printf("%f %f | ", creal(x[i]), cimag(x[i]));
//        printf("%f %f\n", creal(plan->data.filters[0].freq_at(i)), cimag(plan->data.filters[0].freq_at(i)));
//    }
//    printf("%d/%d\n", cnt, n);
//    fftw_free(x);

    int tries = 10;
    int ok = 0;
    for (int iter = 0; iter < tries; ++iter) {

        for (int i = 0; i < N; ++i) {
            out[i][0] = 0;
            out[i][1] = 0;
        }

        for (int i = 0; i < k; ++i) {
            int pos = rand() % (N);
            out[pos][0] = 1;
//            printf("%d ", pos);
        }
//        printf("\n");

        fftw_execute(p);

        for (int i = 0; i < N; ++i) {
            input[i] = in[i][0] + in[i][1] * I;
        }

//        {
//            complex_t* y = (complex_t*) fftw_malloc(n * sizeof(*y));
//            complex_t* v = (complex_t*) fftw_malloc(plan->data.filters[1].B_g * sizeof(*v));
//
//            int sigma = 839;
//            int si = mod_inverse(sigma, n);
//            int b = 250;
//            int a = 0;
//            for (int i = 0; i < n; ++i) {
//                double phi = -((2 * M_PI) / n) * ((sigma * b * i) % n);
//                complex_t value = (cos(phi) + I * sin(phi));
//                y[i] = input[(sigma * (i - a + n)) % n] * value;
//            }
////            for (int i = 0; i < plan->data.filters[1].B_g; ++i) {
////                v[i] = 0;
////                for (int j = 0; j < n / plan->data.filters[1].B_g; ++j) {
////                    v[i] += y[i + plan->data.filters[1].B_g * j];
////                }
////            }
//
////            fftw_dft(v, plan->data.filters[1].B_g, v);
//            fftw_dft(y, n, y);
////            for (int i = 0; i < plan->data.filters[1].B_g; ++i) {
////                printf("%f %f\n", creal(v[i]), cimag(v[i]));
////            }
//            for (int i = 0; i < n; ++i) {
//                if (cabs2(y[i]) > 1e-2) {
//                    printf("%i: %f %f\n", i, creal(y[i]), cimag(y[i]));
//                }
//            }
//            fftw_free(y);
//            fftw_free(v);
//        }

        sfft_output output;

        sfft_exec_multidim(plan, input, &output);
        printf("%ld\n", output.size());
//        for (auto w : output) {
//            printf("%d: %f %f\n", w.first, creal(w.second), cimag(w.second));
//        }
        ok += CheckAnswer(N, output, out);
    }
    printf("%d/%d\n", ok, tries);
    sfft_free(input);
    sfft_free_plan_multidim(plan);

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
    return 0;
}
