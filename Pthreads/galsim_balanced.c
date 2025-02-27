//
// Created by Sam Wallace on 27.02.2025.
//
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

typedef struct {
    int start, end;
    int N;
    int n_threads;
    double *m;
    double *x;
    double *y;
    double eps;
} ThreadData;


double *readData(const char *filename, int N); // Reads the initial condition data [x y m vx vy L]
double *transform(const double *data, int N); // Allocates new 48N bytes in heap for the new order [m][x][y][vx][vy][L]
void SaveLastStep(const char *filename, double *DATA, int N); // Saves final [x y m vx vy L]

void *compute_fores(void *arg) {
    ThreadData *data = (ThreadData *) arg;
    int start = data->start, end = data->end;
    int N = data->N;
    int n_threads = data->n_threads;
    double *m = data->m;
    double *x = data->x;
    double *y = data->y;
    double eps = data->eps;
    double *local_a = calloc(2 * N, sizeof(double));
    double *ax = local_a;
    double *ay = local_a + N;
    double m_j, x_i, y_i;
    int i, j;
    register double F1, F2, dx, dy, dx2, dy2, R2, R, R3, invR3, Gx, Gy, m_i; // Frequently used in innermost loop

    for (i = start; i < end; i += n_threads) {
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
        ax[i] = F1; ay[i] = F2;
    }
    return local_a;
}


int main(int argc, char *argv[]) {

    if (argc != 7) {
        printf("Incorrect usage: ./gaslsim N filname nsteps delta_t graphics n_threads\n");
        return 1;
    }

    int N = atoi(argv[1]);
    const char *filename = argv[2];
    int n_steps = atoi(argv[3]);
    const double dt = atof(argv[4]), G = 100.0 / N, eps = 1e-3;
    // int graphics = atoi(argv[5]);
    int n_threads = atoi(argv[6]);

    // Sanity Check
    printf("n-threads = %d\n", n_threads);

    // Same as Vanilla
    double *data = readData(filename, N);
    double *DATA = transform(data, N);
    double *m = DATA;
    double *x = DATA + N;
    double *y = DATA + 2*N;
    double *u = DATA + 3*N;
    double *v = DATA + 4 * N;
    double *a = calloc(2 * N, sizeof(double));
    double *ax = a;
    double *ay = a + N;
    pthread_t threads[n_threads];
    ThreadData thread_data[n_threads];

    int i, t, n;
    n = 0;


    for (n = 0; n < n_steps; n++) {
        memset(a, 0, 2 * N * sizeof(double));
        for (t = 0; t < n_threads; t++) {
            thread_data[t].start = t; // Start each thread at its index
            thread_data[t].end = N; // End at N particles
            thread_data[t].N = N;
            thread_data[t].m = m;
            thread_data[t].x = x;
            thread_data[t].y = y;
            thread_data[t].eps = eps;
            thread_data[t].n_threads = n_threads;
            pthread_create(&threads[t], NULL, compute_fores, &thread_data[t]);
        }
        // Must wait for all our threads to come back to us for this step iteration
        for (t = 0; t < n_threads; t++) {
            double *local_a = NULL;
            pthread_join(threads[t], (void **) &local_a);
            double *ax_local = local_a;
            double *ay_local = local_a + N;
            for (i = 0; i < N; i++) {
                ax[i] += ax_local[i];
                ay[i] += ay_local[i];
            }
            free(local_a);
        }

        for (i = 0; i < N; i++) {
            v[i] += ax[i] * dt; // Update x-velocity
            x[i] += u[i] * dt; // Update x-position
            v[i] += ay[i] * dt; // Update y-velocity
            y[i] += v[i] * dt; // Update y-position
        }

    }
    SaveLastStep("result_balanced.gal", DATA, N);
    free(a);
    free(data);
    free(DATA);

    // Only update is detroy the mutex -- NO LONGER NECESSSARY with local a's!
    // pthread_mutex_destroy(&mutex);
    return 0;
}

static inline void
Compute_ax_ay(int i, int j, double *__restrict m, double *__restrict x, double *a, double G, double eps) {
    double dx, dy, R2, R, R3, invR3, Gx, Gy;
    // Use particle positions: x[2*i], x[2*i+1] for particle i
    dx = x[2 * j] - x[2 * i];
    dy = x[2 * j + 1] - x[2 * i + 1];

    R2 = dx * dx + dy * dy;
    R = sqrt(R2) + eps;
    R3 = R * R * R;
    invR3 = 1.0 / R3;

    Gx = G * invR3 * dx;
    Gy = G * invR3 * dy;

    // For particle i, add contribution from j
    a[2 * i] += Gx * m[j];
    a[2 * i + 1] += Gy * m[j];
    // For particle j, subtract contribution from i (Newton's third law)
    a[2 * j] -= Gx * m[i];
    a[2 * j + 1] -= Gy * m[i];
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