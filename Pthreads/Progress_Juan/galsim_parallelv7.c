//
// Created by Juan Rodriguez on 2025-03-02.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <pthread.h>
#include "timer.h"

#define eps 1e-3


typedef struct data_structure {
    long rank;
} data_structure;

int N, M, n_steps;
int B_M = 0;
double *m, *x, *y, *u, *v, *L, *a, *ax, *ay, *a_loc;
pthread_mutex_t b_mutex;
pthread_cond_t b_cond_var;
double G = 0, dt = 1e-5;

double *readData(const char *filename, int N);
double *transform(const double *data, int N);
void SaveLastStep(const char *filename, double *DATA, int N);
void printArray1D(double *v, int N);
void *thread_function(void *rank);
void Barrier_init(void);
void Barrier(void);
void Barrier_destroy(void);

int main(int argc, char *argv[]) {
    N = atoi(argv[1]);
    const char *filename = argv[2];
    n_steps = atoi(argv[3]);
    dt = atof(argv[4]);
    G = 100.0 / N;
    M = atoi(argv[6]);

    double *data = readData(filename, N), *DATA = transform(data, N);
    m = DATA, x = DATA + N, y = DATA + 2 * N;
    u = DATA + 3 * N, v = DATA + 4 * N, L = DATA + 5* N;
    a = calloc(2 * N, sizeof(double));
    a_loc = calloc(M * 2 * N,  sizeof(double));
    ax = a;
    ay = a+N;
    double start, end;
    long t;
    pthread_t threads[M];
    Barrier_init();
    GET_TIME(start);
    for (t = 0; t < M; t++)
        pthread_create(&threads[t], NULL, thread_function, (void *) t);
    for (t = 0; t < M; t++)
        pthread_join(threads[t], NULL);
    GET_TIME(end);
    printf("Time elapsed: %lf\n", end - start);
    SaveLastStep("../compare_gal_files/result_balanced.gal", DATA, N);
    free(a);
    free(data);
    free(DATA);
}

void *thread_function(void *rank) {
    int rank_ = (int)((long) rank);
    int t, n, i;
    int bfirst, blast, bincr, cfirst, clast, cincr, iterations;

    cfirst = rank_, clast = N, cincr = M;

    int work_ratio = N / M;
    int remainder = N % M;
    bincr = 1;
    if (rank_ < remainder){
        iterations = work_ratio + 1;
        bfirst = rank_ * iterations;
    } else {
        iterations = work_ratio;
        bfirst = rank_ * iterations + remainder;
    }
    blast = bfirst + iterations;

    double *ax_loc = a_loc + rank_ * 2 * N;
    double *ay_loc = ax_loc + N;
    for (n = 0; n < n_steps; n++) {
        memset(ax_loc, 0, N * sizeof(double));
        memset(ay_loc, 0, N * sizeof(double));

        double m_j, x_i, y_i;
        register double F1, F2, dx, dy, dx2, dy2, R2, R, R3, invR3, Gx, Gy, m_i;
        Barrier();
        for (i = cfirst; i < clast; i += cincr) {
            x_i = x[i];
            y_i = y[i];
            m_i = m[i];
            F1 = ax_loc[i];
            F2 = ay_loc[i];
#pragma GCC ivdep
            for (int j = i + 1; j < N; j++) {
                dx = x[j] - x_i;
                dy = y[j] - y_i;
                dx2 = dx * dx;
                dy2 = dy * dy;
                R2 = dx2 + dy2;
                R = sqrt(R2) + eps;
                R3 = R * R * R;
                invR3 = 1.0 / R3;
                Gx = invR3 * dx;
                Gy = invR3 * dy;
                ax_loc[j] -= Gx * m_i;
                ay_loc[j] -= Gy * m_i;
                m_j = m[j];
                F1 += Gx * m_j;
                F2 += Gy * m_j;
            }
            ax_loc[i] = F1;
            ay_loc[i] = F2;
        }
        Barrier();
        for (i = bfirst; i < blast; i += bincr) {
            //ax[i] = ay[i] = 0.0;
            register double F1, F2;
            F1 = F2 = 0.0;
            int idx = i;
#pragma GCC ivdep
            for (t = 0; t < M; t++){
                //ax[idx] += a_loc[t*(2*N) + idx];
                //ay[idx] += a_loc[t*(2*N) + N + idx];
                F1 += a_loc[t*(2*N) + idx];
                F2 += a_loc[t*(2*N) + N + idx];
            }
            u[i] += G * (F1 * dt);
            v[i] += G * (F2 * dt);
            x[i] += u[i] * dt;
            y[i] += v[i] * dt;
        }
        //Barrier();
//        for (i = bfirst; i < blast; i += bincr){
//            u[i] += G * (ax[i] * dt);
//            x[i] += u[i] * dt;
//            v[i] += G * (ay[i] * dt);
//            y[i] += v[i] * dt;
//        }
        Barrier();
    }
    return NULL;
}

void Barrier_init(void) {
    B_M = 0;
    pthread_mutex_init(&b_mutex, NULL);
    pthread_cond_init(&b_cond_var, NULL);
}

void Barrier(void) {
    pthread_mutex_lock(&b_mutex);
    B_M++;
    if (B_M == M) {
        B_M = 0;
        pthread_cond_broadcast(&b_cond_var);
    } else {
        while (pthread_cond_wait(&b_cond_var, &b_mutex) != 0);
    }
    pthread_mutex_unlock(&b_mutex);
}

void Barrier_destroy(void) {
    pthread_mutex_destroy(&b_mutex);
    pthread_cond_destroy(&b_cond_var);
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

void printArray1D(double *v, int N) {
    for (int i = 0; i < N; i++) {
        printf("%lf\n", v[i]);
    }
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
    double *L = DATA + 5 * N;
    for (int i = 0; i < N; i++) {
        fwrite(&x[i], sizeof(double), 1, file);
        fwrite(&y[i], sizeof(double), 1, file);
        fwrite(&m[i], sizeof(double), 1, file);
        fwrite(&u[i], sizeof(double), 1, file);
        fwrite(&v[i], sizeof(double), 1, file);
        fwrite(&L[i], sizeof(double), 1, file);
    }
    fclose(file);
}
