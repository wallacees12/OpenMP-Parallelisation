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

typedef struct thread_container {
    int start, end;
    double *m, *x, *y, *u, *v;
} thread_container;

typedef struct data_structure {
    double *m, *x, *y, *u, *v;
} data_structure;

double G = 0, dt = 1e-5;
int N, M, n_steps;

double *readData(const char *filename, int N); // Reads the initial condition data [x y m vx vy L]
double *transform(const double *data, int N); // Allocates new 48N bytes in heap for the new order [m][x][y][vx][vy][L]
void SaveLastStep(const char *filename, double *DATA, int N); // Saves final [x y m vx vy L]
void printArray1D(double *v, int N);

void *compute_accelerations(void *arg);

static inline void InitializeContainers(thread_container *thread_containers, data_structure *data);

int main(int argc, char *argv[]) {

    N = atoi(argv[1]);
    const char *filename = argv[2];
    n_steps = atoi(argv[3]);
    dt = atof(argv[4]);
    G = 100.0 / N;
    M = atoi(argv[6]);

    double *data = readData(filename, N), *DATA = transform(data, N);

    double *m = DATA, *x = DATA + N, *y = DATA + 2 * N;
    double *u = DATA + 3 * N, *v = DATA + 4 * N;
    double *a = calloc(2 * N, sizeof(double));
    double *ax = a, *ay = a + N;
    double start, end;

    data_structure dStructure = {m, x, y, u, v};
    pthread_t threads[M];
    thread_container thread_containers[M];
    InitializeContainers(thread_containers, &dStructure);

    int n;
    GET_TIME(start);
    for (n = 0; n < n_steps; n++) {
        // Iteratively create and join threads -> Increases thread creation overhead ->
        // -> For larger n_steps -> Affects performance
        for (int i = 0; i < M; i++)
            pthread_create(&threads[i], NULL, compute_accelerations, &thread_containers[i]);

        for (int t = 0; t < M; t++) {
            double *local_a;
            pthread_join(threads[t], (void **) &local_a);
            double *ax_local = local_a, *ay_local = local_a + N;
            for (int i = 0; i < N; i++) {
                ax[i] += ax_local[i];
                ay[i] += ay_local[i];
            }
            free(local_a);
        }
        // This section could further be improved for better parallelism
        // N is comparatively small to N(N-1)/2 for larger N
        for (int i = 0; i < N; i++) {
            u[i] += G * (ax[i] * dt);
            x[i] += u[i] * dt;
            v[i] += G * (ay[i] * dt);
            y[i] += v[i] * dt;
            ax[i] = ay[i] = 0.0; // To avoid memset. Turns out to be faster maybe?
        }
    }
    GET_TIME(end);
    printf("Time ellapsed: %lf\n", end - start);

    SaveLastStep("../compare_gal_files/result_balanced.gal", DATA, N);
    free(a);
    free(data);
    free(DATA);
}

void *compute_accelerations(void *arg) {
    thread_container *thread_data = (thread_container *) arg;
    int start = thread_data->start;
    int end = thread_data->end;

    double *m = thread_data->m, *x = thread_data->x, *y = thread_data->y;
    double *a = calloc(2 * N, sizeof(double));
    double *ax = a, *ay = a + N;
    int i, j;
    double m_j, x_i, y_i;
    register double F1, F2, dx, dy, dx2, dy2, R2, R, R3, invR3, Gx, Gy, m_i; // Frequently used in innermost loop
    for (i = start; i < end; i += M) {
        x_i = x[i];
        y_i = y[i];
        m_i = m[i];
        F1 = ax[i];
        F2 = ay[i];
#pragma GCC ivdep
        for (j = i + 1; j < N; j++) {
            dx = x[j] - x_i;
            dy = y[j] - y_i;
            dx2 = dx * dx;
            dy2 = dy * dy;
            R2 = dx2 + dy2;
            R = sqrt(R2) + eps;
            R3 = R * R;
            R3 *= R;
            invR3 = 1.0 / R3;
            Gx = Gy = invR3;
            Gx *= dx;
            Gy *= dy;
            ax[j] -= Gx * m_i;
            ay[j] -= Gy * m_i;
            m_j = m[j];
            F1 += Gx * m_j;
            F2 += Gy * m_j;
        }
        ax[i] = F1;
        ay[i] = F2;
    }
    return a;
}


static inline void InitializeContainers(thread_container *thread_containers, data_structure *data) {
    double *m = data->m, *x = data->x, *y = data->y;
    double *u = data->u, *v = data->v;
    for (int i = 0; i < M; i++) {
        thread_containers[i].start = i;
        thread_containers[i].end = N;
        thread_containers[i].m = m;
        thread_containers[i].x = x;
        thread_containers[i].y = y;
        thread_containers[i].u = u;
        thread_containers[i].v = v;
    }
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
        DATA[i] = data[6 * i + 2]; // Since this function is called only once defining idx = 6*i has little effect
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
    // Define pointers to segments in DATA for more readable write instructions
    double *m = DATA;
    double *x = DATA + N;
    double *y = DATA + 2 * N;
    double *u = DATA + 3 * N;
    double *v = DATA + 4 * N;
    double *L = DATA + 5 * N;

    for (int i = 0; i < N; i++) {
        fwrite(&x[i], sizeof(double), 1, file); // Write 1 double for each quantity
        fwrite(&y[i], sizeof(double), 1, file);
        fwrite(&m[i], sizeof(double), 1, file);
        fwrite(&u[i], sizeof(double), 1, file);
        fwrite(&v[i], sizeof(double), 1, file);
        fwrite(&L[i], sizeof(double), 1, file);
    }
    fclose(file);
}