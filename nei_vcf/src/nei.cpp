#ifdef _WIN32
#define LIBRARY_API extern "C" __declspec(dllexport)
#else
#define LIBRARY_API extern "C"
#endif

#include<math.h>

static inline double nei_dist_single(double* a, double* b) {
    return 1 - (
        sqrt(a[0] * b[0]) +
        sqrt(a[1] * b[1]) +
        sqrt(a[2] * b[2]) +
        sqrt(a[3] * b[3])
    );
}

static inline double sum_4(double *x) {
    return x[0] + x[1] + x[2] + x[3];
}

LIBRARY_API void nei_dist(double* variants, int* shape, double* distances_summed, int* counts) {
    for (int k = 0; k < shape[0]; k++) { // k (variants) first to avoid cache misses
        for (int i = 0; i < shape[1]; i++) {
            for (int j = i + 1; j < shape[1]; j++) {
            
                // printf("i: %d; j:%d; k: %d\n", i, j, k);
                double *a = variants + (i * 4 + k * shape[1] * 4);
                double *b = variants + (j * 4 + k * shape[1] * 4);

                if (sum_4(a) == 0 || sum_4(b) == 0) continue;

                double dist = nei_dist_single(a, b);
                distances_summed[j + i * shape[1]] += dist;
                distances_summed[i + j * shape[1]] += dist;
                counts[j + i * shape[1]] += 1;
                counts[i + j * shape[1]] += 1;
            }
        }
    }
    for (int i = 0; i < shape[1]; i++) {
        counts[i + i * shape[1]] = shape[0];
    }
}
