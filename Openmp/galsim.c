#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>

typedef struct{
    int N;
    double *m;
    double *x;
    double G;
    double eps;
} ThreadedData;

static inline void Compute_ax_ay(int i, int j, double *__restrict m,
                                 double *__restrict x, double *ax_i, double *ay_i,
                                 double *ax_j, double *ay_j, double G, double eps);

double *readData(const char *filename, int N);
double *transform(const double *data, int N);
void SaveLastStep(const char *filename, double *DATA, int N);

void compute_forces(ThreadedData *data, double *a, int n_threads) {
    int N = data->N;
    double *m = data->m;
    double *x = data->x;
    double G = data->G;
    double eps = data->eps;

    // Create a per-thread local acceleration array
    double **local_a = (double **) malloc(n_threads * sizeof(double *));
    #pragma omp parallel num_threads(n_threads)
    {
        int tid = omp_get_thread_num();
        local_a[tid] = (double *) calloc(2 * N, sizeof(double));

        #pragma omp for schedule(dynamic)
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                double ax_i = 0, ay_i = 0, ax_j = 0, ay_j = 0;
                Compute_ax_ay(i, j, m, x, &ax_i, &ay_i, &ax_j, &ay_j, G, eps);

                local_a[tid][2 * i]     += ax_i;
                local_a[tid][2 * i + 1] += ay_i;
                local_a[tid][2 * j]     += ax_j;
                local_a[tid][2 * j + 1] += ay_j;
            }
        }

        // Merge local results back into global `a[]`
        #pragma omp critical
        for (int i = 0; i < 2 * N; i++) {
            a[i] += local_a[tid][i];
        }

        free(local_a[tid]);
    }

    free(local_a);
}

int main(int argc, char *argv[]) {
    if (argc != 7) {
        printf("Incorrect usage: ./gaslsim N filename nsteps delta_t graphics n_threads\n");
        return 1;
    }

    int N = atoi(argv[1]);
    const char *filename = argv[2];
    int n_steps = atoi(argv[3]);
    double dt = atof(argv[4]), G = 100.0 / N, eps = 1e-3;
    int n_threads = atoi(argv[6]);

    double *data = readData(filename, N);
    double *DATA = transform(data, N);
    double *a = calloc(2 * N, sizeof(double));
    double *m = DATA;
    double *x = DATA + N;
    double *v = DATA + 3 * N;

    ThreadedData threadedData = {N, m, x, G, eps};
    for (int n = 0; n < n_steps; n++) {
        memset(a, 0, 2 * N * sizeof(double));
        compute_forces(&threadedData, a, n_threads);

        #pragma omp parallel for num_threads(n_threads)
        for (int i = 0; i < 2 * N; i += 2) {
            v[i] += a[i] * dt;
            x[i] += v[i] * dt;
            v[i + 1] += a[i + 1] * dt;
            x[i + 1] += v[i + 1] * dt;
        }
    }
    SaveLastStep("result.gal", DATA, N);
    free(a);
    free(data);
    free(DATA);
    return 0;
}

static inline void Compute_ax_ay(int i, int j, double *__restrict m,
                                 double *__restrict x, double *ax_i, double *ay_i,
                                 double *ax_j, double *ay_j, double G, double eps) {
    double dx = x[2 * j] - x[2 * i];
    double dy = x[2 * j + 1] - x[2 * i + 1];

    double R2 = dx * dx + dy * dy;
    double R = sqrt(R2) + eps;
    double R3 = R * R * R;
    double invR3 = 1.0 / R3;

    double Gx = G * invR3 * dx;
    double Gy = G * invR3 * dy;

    *ax_i = Gx * m[j];
    *ay_i = Gy * m[j];
    *ax_j = -Gx * m[i];
    *ay_j = -Gy * m[i];
}

double *readData(const char *filename, int N) {
    FILE *file = fopen(filename, "rb");

    if (file == NULL) {
        printf("Error opening file %s\n", filename);
        exit(1);
    }

    double *data = malloc(6 * N * sizeof(double));

    if (fread(data, sizeof(double), 6 * N, file) != 6 * N) {
        printf("Error reading file %s\n", filename);
        fclose(file);
        free(data);
        exit(1);
    }

    fclose(file);
    return data;
}


double *transform(const double *data, int N) {
    double *DATA = malloc(6 * N * sizeof(double));
    for (int i = 0; i < N; i++) {
        DATA[i] = data[6 * i + 2]; // m0, m1, ..., m(N-1),
        DATA[N + 2 * i] = data[6 * i]; // m(N-1), x0,
        DATA[N + 2 * i + 1] = data[6 * i + 1]; // m(N-1), x0, y0
        DATA[3 * N + 2 * i] = data[6 * i + 3];
        DATA[3 * N + 2 * i + 1] = data[6 * i + 4];
        DATA[5 * N + i] = data[6 * i + 5];
    }
    return DATA;
}


void SaveLastStep(const char *filename, double *DATA, int N) {
    FILE *file = fopen(filename, "wb");
    if (file == NULL) {
        printf("Error opening file %s\n", filename);
        exit(1);
    }
    double *m = DATA;
    double *x = DATA + N;
    double *v = DATA + 3 * N;
    double *L = DATA + 5 * N;
    for (int i = 0; i < N; i++) {
        fwrite(&x[i * 2], sizeof(double), 1, file); // &x[] pointer to accessed value
        fwrite(&x[i * 2 + 1], sizeof(double), 1, file);
        fwrite(&m[i], sizeof(double), 1, file);
        fwrite(&v[i * 2], sizeof(double), 1, file);
        fwrite(&v[i * 2 + 1], sizeof(double), 1, file);
        fwrite(&L[i], sizeof(double), 1, file);
    }
    fclose(file);
}