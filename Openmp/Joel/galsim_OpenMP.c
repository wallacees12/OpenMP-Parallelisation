#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>


double *readData(const char *filename, int N);

double *transform(const double *data, int N);

void alloc_a(double*** a, int nThreads, int N){
    int i;
    *a = (double**)malloc(nThreads*sizeof(double));
    for(i=0;i<nThreads;i++){
        (*a)[i] = (double*)calloc(N, sizeof(double));
    }
}

void free_a(double** a, int nThreads){
    int i;
    for(i = 0; i<nThreads; i++){
        free(a[i]);
    }
    free(a);
}

void SaveLastStep(const char *filename, double *DATA, int N); // Doesn't step out of a mountain

int main(int argc, char *argv[]) {

    if (argc < 7) {
        printf("Too few of arguments.\n");
        return 1;
    } else if (7 < argc) {
        printf("Too many arguments.\n");
        return 1;
    }

    int N = atoi(argv[1]);
    const char *filename = argv[2];
    int n_steps = atoi(argv[3]);
    const double dt = atof(argv[4]), G = 100.0 / N, eps = 1e-3;
    int graphics = atoi(argv[5]);
    int nThreads = atoi(argv[6]);

    double *data = readData(filename, N);
    double *DATA = transform(data, N);
    // [m][x][y][u][v] ~ [ax][ay]
    double *a = calloc(2 * N, sizeof(double));
    double **local_ax;
    double **local_ay;
    double *m = DATA;
    double *x = DATA + N;
    double *y = DATA + 2 * N;
    double *u = DATA + 3 * N;
    double *v = DATA + 4 * N;
    double *ax = a;
    double *ay = a + N;

    double m_j, x_i, y_i;
    int i, j, n;
    register double F1, F2, dx, dy, dx2, dy2, R2, R, R3, invR3, Gx, Gy, m_i;
    n = 0;
    alloc_a(&local_ax, nThreads, N);
    alloc_a(&local_ay, nThreads, N);

    while (n < n_steps) {
        #pragma omp paralllel num_threads(nThreads)
        {
            int thrid = omp_get_thread_num();

            #pragma omp parallel for schedule
            for (i = 0; i < N; i++) {
                x_i = x[i]; y_i = y[i]; m_i = m[i];
                F1 = ax[i]; F2 = ay[i];
                #pragma GCC ivdep
                for (j = i + 1; j < N; j++) {
                    dx = x[j] - x_i; dy = y[j] - y_i;
                    dx2 = dx * dx; dy2 = dy * dy;
                    R2 = dx2 + dy2;
                    R = sqrt(R2) + eps;
                    R3 = R * R;
                    R3 *= R;
                    invR3 = 1.0 / R3;
                    Gx = Gy = invR3;
                    Gx *= dx; Gy *= dy;
                    ax[j] -= Gx * m_i; ay[j] -= Gy * m_i;
                    m_j = m[j];
                    F1 += Gx * m_j; F2 += Gy * m_j;
                }
                u[i] += G * (F1 * dt);
                x[i] += u[i] * dt;
                v[i] += G * (F2 * dt);
                y[i] += v[i] * dt;
            }
            n++;
            memset(a, 0, 2 * N * sizeof(double));
        }
    }
    SaveLastStep("../compare_gal_files/result.gal", DATA, N);
    free(a);
    free_a(local_ax,nThreads);
    free_a(local_ay,nThreads);
    free(data);
    free(DATA);
    return 0;
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
        DATA[i] = data[6 * i + 2];
        DATA[N + i] = data[6 * i];
        DATA[2 * N + i] = data[6 * i + 1];
        DATA[3 * N + i] = data[6 * i + 3];
        DATA[4 * N + i] = data[6 * i + 4];
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
    double *y = DATA + 2 * N;
    double *u = DATA + 3 * N;
    double *v = DATA + 4 * N;
    double *L = DATA + 5 * N;

    for (int i = 0; i < N; i++) {
        fwrite(&x[i], sizeof(double), 1, file); // &x[] pointer to accessed value
        fwrite(&y[i], sizeof(double), 1, file);
        fwrite(&m[i], sizeof(double), 1, file);
        fwrite(&u[i], sizeof(double), 1, file);
        fwrite(&v[i], sizeof(double), 1, file);
        fwrite(&L[i], sizeof(double), 1, file);
    }
    fclose(file);
}