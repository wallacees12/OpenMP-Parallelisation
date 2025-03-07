//
// Created by Juan Rodriguez on 2025-03-03.
//
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include "timer.h"


int n_threads = 8;

void transform(double *data, int N) {
    double org_data[6 * N];

    for (int i = 0; i < N; i++) {
        org_data[i] = data[6 * i + 2];

        org_data[N + i] = data[6 * i];
        org_data[2 * N + i] = data[6 * i + 1];

        org_data[3 * N + i] = data[6 * i + 3];
        org_data[4 * N + i] = data[6 * i + 4];

        org_data[5 * N + i] = data[6 * i + 5];
    }

    memmove(data, org_data, sizeof(org_data));
}

void Detransform(double *org_data, int N) {
    double data[6 * N];
    for (int i = 0; i < N; i++) {
        data[6 * i + 2] = org_data[i];

        data[6 * i] = org_data[N + i];
        data[6 * i + 1] = org_data[2 * N + i];

        data[6 * i + 3] = org_data[3 * N + i];
        data[6 * i + 4] = org_data[4 * N + i];

        data[6 * i + 5] = org_data[5 * N + i];
    }
    memmove(org_data, data, sizeof(data));
}

const float circleRadius = 0.002f, circleColor = 0;
const int windowWidth = 800;

//--------------------------


pthread_mutex_t my_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t my_cond = PTHREAD_COND_INITIALIZER;
int th_counter = 0;


typedef struct thread {
    int n_steps;
    int N;
    double Dt;

    int id;
    int start;
    int jump;
    int end;
    int idStart;
    int idEnd;

    double *m;
    double *x;
    double *a;
    double *v;

    double *a_threads;
} thread;


void *threadfun(void *localThread) {

    thread *th = (thread *) localThread;

    int n_steps = th->n_steps;
    int N = th->N;
    double Dt = th->Dt;
    int id = th->id;
    int start = th->start;
    int jump = th->jump;
    int end = th->end;
    int idStart = th->idStart;
    int idEnd = th->idEnd;

    double *m = th->m;
    double *x = th->x;
    double *a = th->a;
    double *vx = th->v;
    double *vy = vx + N;


    double *a_threads = (th->a_threads); //Address of contiguously allocated a_threads for each thread.
    double *a_th = (th->a_threads) + 2 * N * id;

    double *ax_th = a_th;
    double *ay_th = a_th + N;
    double *y = x + N;

    register double Dx, Dx2, Dy, Dy2, gx, gy, R, r2, Rp3, inv_R3, Fxi, Fyi, xi, yi, mi;
    const double G = (double) 100.0 / N;
    const double eps = 1e-3;

    register int n = 0;
    while (n < n_steps) {
        memset(a_th, 0, 2 * N * sizeof(double));
        for (int i = start; i < end; i += jump) {
            Fxi = ax_th[i];
            Fyi = ay_th[i];
            xi = x[i];
            yi = y[i];
            mi = m[i];
#pragma GCC ivdep
            for (int j = i + 1; j < N; j++) {
                Dx = x[j] - xi;
                Dy = y[j] - yi;
                Dx2 = Dx * Dx;
                Dy2 = Dy * Dy;
                r2 = Dx2 + Dy2;
                R = sqrt(r2) + eps;
                Rp3 = R * R * R;
                inv_R3 = 1 / (Rp3);

                gx = Dx * inv_R3;
                gy = Dy * inv_R3;

                // a_x[i], a_x[j]

                Fxi += gx * m[j];
                ax_th[j] += -gx * mi;

                // a_y[i], a_y[j]

                Fyi += gy * m[j];
                ay_th[j] += -gy * mi;
            }

            ax_th[i] = Fxi;
            ay_th[i] = Fyi;
        }

        pthread_mutex_lock(&my_mutex);
        th_counter += 1;
        if (th_counter < n_threads) {
            pthread_cond_wait(&my_cond, &my_mutex);
        }
        if (th_counter == n_threads) {
            th_counter = 0;
            pthread_cond_broadcast(&my_cond);
        }
        pthread_mutex_unlock(&my_mutex);


#pragma GCC ivdep
        for (int i = idStart; i < idEnd; i++) {
            register double Fxk = 0;
            register double Fyk = 0;
            for (int k = 0; k < n_threads; k++) { // Loop unrolling maybe. Cache access concern: Probably fine though
                Fxk += a_threads[k * 2 * N + i];
                Fyk += a_threads[N * (2 * k + 1) + i];
            }
            vx[i] += G * (Fxk * Dt);
            vy[i] += G * (Fyk * Dt);
            x[i] += vx[i] * Dt;
            y[i] += vy[i] * Dt;
        }

        n++;
        pthread_mutex_lock(&my_mutex);
        th_counter += 1;
        if (th_counter < n_threads) {
            pthread_cond_wait(&my_cond, &my_mutex);
        } else {
            th_counter = 0;
            pthread_cond_broadcast(&my_cond);
        }
        pthread_mutex_unlock(&my_mutex);
    }
    return NULL;
}


int main(int argc, char *argv[]) {

    int N = atoi(argv[1]);
    char *file_name = argv[2];
    char *result_name = "../compare_gal_files/result_parallel.gal";// result.gal
    int n_steps = atoi(argv[3]);

    double Dt = strtod(argv[4], NULL);
    int graphics = atoi(argv[5]);
    n_threads = atoi(argv[6]);
    const double G = (double) 100.0 / N;
    const double eps = 1e-3;
    float L = 1, W = 1;

    printf("Press 'q' to quit\n");

    FILE *fp = fopen(file_name, "rb");

    //De momento guardo memoria en el stack.
    double *data = malloc(6 * N * sizeof(double));

    // Open file
    fread(data, sizeof(double), 6 * N, fp);

    //Reorganize Data
    transform(data, N);

    //Close file
    fclose(fp);

    double *a = malloc(2 * N * sizeof(double));
    double *ax = a;
    double *ay = a + N;

    double *m = data;
    double *x = data + N;
    double *y = data + 2 * N;
    double *vx = data + 3 * N;
    double *vy = data + 4 * N;


    const int chunk = N / n_threads;
    thread mythread[n_threads];
    double *a_threads = (double *) malloc(n_threads * 2 * N * sizeof(double));
    for (int k = 0; k < n_threads; k++) {
        mythread[k].n_steps = n_steps;
        mythread[k].N = N;
        mythread[k].Dt = Dt;

        mythread[k].id = k;
        mythread[k].start = k;
        mythread[k].jump = n_threads;
        mythread[k].end = N;

        mythread[k].idStart = k * chunk;
        mythread[k].idEnd = (k + 1) * chunk;
        if (k == n_threads - 1) {
            mythread[k].idEnd = N;
        }

        mythread[k].m = m;
        mythread[k].x = x;
        mythread[k].a = a;
        mythread[k].v = vx;
        mythread[k].a_threads = a_threads;
    }

    pthread_t threads[n_threads];
    double start, end;
    GET_TIME(start);
    for (int k = 0; k < n_threads; k++) {
        pthread_create(&(threads[k]), NULL, threadfun, &(mythread[k]));
    }

    for (int k = 0; k < n_threads; k++) {
        pthread_join(threads[k], NULL);
    }
    GET_TIME(end);
    printf("Time elapsed: %lf\n", end - start);

    Detransform(data, N);

    FILE *file_write;
    file_write = fopen(result_name, "wb");
    if (fp == NULL) {
        printf("Error opening file %s\n", result_name);
        exit(1);
    }
    fwrite(data, sizeof(double), 6 * N, file_write);

    fclose(file_write);

    free(data);
    free(a_threads);
    free(a);

    return 0;
}
// gcc -Wall -o galsim_parallelv8 galsim_parallelv8.c -lm -lpthread
// ./galsim_parallelv8 3000 ./input_data/ellipse_N_03000.gal 100 1e-5 0 4