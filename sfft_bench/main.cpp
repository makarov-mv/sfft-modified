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

    sfft_plan_multidim* plan = sfft_make_plan_multidim(n, 1, k, 2);
    complex_t* sin = (complex_t*) sfft_malloc(sizeof(complex_t) * n);

    int tries = 1;
    int ok = 0;
    for (int iter = 0; iter < tries; ++iter) {

        for (int i = 0; i < n; ++i) {
            out[i][0] = 0;
            out[i][1] = 0;
        }

        for (int i = 0; i < k; ++i) {
            int pos = rand() % n;
            out[pos][0] = 1;
//            printf("%d ", pos);
        }
//        printf("\n");

        fftw_execute(p);



        for (int i = 0; i < n; ++i) {
            sin[i] = in[i][0] + in[i][1] * I;
        }


        sfft_output output;

        sfft_exec_multidim(plan, sin, &output);
//        printf("%ld\n", output.size());
//        for (auto w : output) {
//            printf("%d: %f %f\n", w.first, creal(w.second), cimag(w.second));
//        }
        ok += CheckAnswer(n, output, out);
    }
    printf("%d/%d\n", ok, tries);
    sfft_free(sin);
//    sfft_free_plan(plan);

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
    return 0;
}