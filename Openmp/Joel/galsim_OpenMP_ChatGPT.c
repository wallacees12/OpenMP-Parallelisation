#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>

double *readData(const char *filename, int N);
double *transform(const double *data, int N);
void SaveLastStep(const char *filename, double *DATA, int N);

void alloc_a(double ***a, int nThreads, int N) {
    *a = (double **) malloc(nThreads * sizeof(double *));
    if (!*a) {
        printf("Memory allocation failed for `a`\n");
        exit(1);
    }
    for (int i = 0; i < nThreads; i++) {
        (*a)[i] = (double *) calloc(N, sizeof(double));
        if (!(*a)[i]) {
            printf("Memory allocation failed for `a[%d]`\n", i);
            exit(1);
        }
    }
}

void free_a(double **a, int nThreads) {
    for (int i = 0; i < nThreads; i++) {
        free(a[i]);
    }
    free(a);
}

int main(int argc, char *argv[]) {
    if (argc != 7) { 
        printf("Usage: %s N filename n_steps dt graphics nThreads\n", argv[0]);
        return 1;
    }

    int N = atoi(argv[1]);
    const char *filename = argv[2];
    int n_steps = atoi(argv[3]);
    double dt = atof(argv[4]), G = 100.0 / N, eps = 1e-3;
    int graphics = atoi(argv[5]);
    int nThreads = atoi(argv[6]);

    omp_set_num_threads(nThreads);  

    double *data = readData(filename, N);
    if (!data) {
        printf("Error: Unable to read input data.\n");
        return 1;
    }

    double *DATA = transform(data, N);
    if (!DATA) {
        free(data);
        printf("Error: Data transformation failed.\n");
        return 1;
    }

    double *a = (double *) calloc(2 * N, sizeof(double));
    if (!a) {
        printf("Memory allocation failed for `a`.\n");
        free(data);
        free(DATA);
        return 1;
    }

    double **local_ax, **local_ay;
    alloc_a(&local_ax, nThreads, N);
    alloc_a(&local_ay, nThreads, N);

    double *m = DATA, *x = DATA + N, *y = DATA + 2 * N;
    double *u = DATA + 3 * N, *v = DATA + 4 * N;
    double *ax = a, *ay = a + N;

    int n = 0;
    while (n < n_steps) {
        
        /* Reset forces */
        memset(ax, 0, N * sizeof(double));
        memset(ay, 0, N * sizeof(double));

        /* Compute forces */
#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < N; i++) {
            int thrid = omp_get_thread_num();
            double x_i = x[i], y_i = y[i], m_i = m[i];
            double A1 = 0.0, A2 = 0.0;

            for (int j = 0; j < N; j++) {
                if (i == j) continue;

                double dx = x[j] - x_i;
                double dy = y[j] - y_i;
                double R2 = dx * dx + dy * dy;
                double R = sqrt(R2) + eps;
                double invR3 = 1.0 / (R * R * R);
                double Gx = invR3 * dx, Gy = invR3 * dy;

                A1 += Gx * m[j];
                A2 += Gy * m[j];
            }
            local_ax[thrid][i] = A1;
            local_ay[thrid][i] = A2;
        }

        /* Accumulate results safely */
#pragma omp parallel for
        for (int i = 0; i < N; i++) {
            for (int t = 0; t < nThreads; t++) {
                ax[i] += local_ax[t][i];
                ay[i] += local_ay[t][i];
            }
        }

        /* Update positions */
#pragma omp parallel for
        for (int i = 0; i < N; i++) {
            u[i] += G * ax[i] * dt;
            x[i] += u[i] * dt;
            v[i] += G * ay[i] * dt;
            y[i] += v[i] * dt;
        }

        /* Reset per-thread storage */
#pragma omp parallel
        {
            int thrid = omp_get_thread_num();
            memset(local_ax[thrid], 0, N * sizeof(double));
            memset(local_ay[thrid], 0, N * sizeof(double));
        }
        
        n++;
    }

    SaveLastStep("galsim_OpenMP.gal", DATA, N);

    free(a);
    free_a(local_ax, nThreads);
    free_a(local_ay, nThreads);
    free(data);
    free(DATA);

    return 0;
}

/* Safe File Reading */
double *readData(const char *filename, int N) {
    FILE *file = fopen(filename, "rb");
    if (!file) {
        printf("Error opening file: %s\n", filename);
        return NULL;
    }

    double *data = (double *) malloc(6 * N * sizeof(double));
    if (!data) {
        printf("Memory allocation failed for data.\n");
        fclose(file);
        return NULL;
    }

    if (fread(data, sizeof(double), 6 * N, file) != 6 * N) {
        printf("Error reading file: %s\n", filename);
        fclose(file);
        free(data);
        return NULL;
    }

    fclose(file);
    return data;
}

/* Transform Data */
double *transform(const double *data, int N) {
    double *DATA = (double *) malloc(6 * N * sizeof(double));
    if (!DATA) {
        printf("Memory allocation failed for transformed data.\n");
        return NULL;
    }

    for (int i = 0; i < N; i++) {
        DATA[i] = data[6 * i + 2];
        DATA[N + i] = data[6 * i];
        DATA[2 * N + i] = data[6 * i + 1];
        DATA[3 * N + i] = data[6 * i + 3];
        DATA[4 * N + i] = data[6 * i + 4];
        DATA[5 * N + i] = data[6 * i + 5];
    }
    return DATA;
}

/* Save Final State */
void SaveLastStep(const char *filename, double *DATA, int N) {
    FILE *file = fopen(filename, "wb");
    if (!file) {
        printf("Error opening file: %s\n", filename);
        exit(1);
    }

    for (int i = 0; i < N; i++) {
        fwrite(&DATA[N + i], sizeof(double), 1, file); // x
        fwrite(&DATA[2 * N + i], sizeof(double), 1, file); // y
        fwrite(&DATA[i], sizeof(double), 1, file); // m
        fwrite(&DATA[3 * N + i], sizeof(double), 1, file); // u
        fwrite(&DATA[4 * N + i], sizeof(double), 1, file); // v
        fwrite(&DATA[5 * N + i], sizeof(double), 1, file); // L
    }

    fclose(file);
}