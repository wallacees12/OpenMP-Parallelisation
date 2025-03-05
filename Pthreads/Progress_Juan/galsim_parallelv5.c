//
// Created by Juan Rodriguez on 2025-03-02.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include "timer.h"

#define X 0
#define Y 1
#define eps 1e-3
#define DIM 2

double G = 0, dt = 1e-5;
const int BLOCK = 0;
const int CYCLIC = 1;
typedef double vect_t[DIM];

typedef struct data_structure {
    double *m, *x, *y, *u, *v;
} data_structure;

struct particle_s {
    double m;
    vect_t position;
    vect_t velocity;
};

int M;
int N;
int n_steps;
struct particle_s *current;
vect_t *a;
vect_t *loc_a;
int B_M = 0;
pthread_mutex_t b_mutex;
pthread_cond_t b_cond_var;

void Loop_schedule(int my_rank, int M, int N, int sched, int *first_p, int *last_pm, int *incr_p);

void *thread_function(void *rank);

void Compute_force(int part, vect_t forces[]);

void Update_part(int part);

void Barrier_init(void);

void Barrier(void);

void Barrier_destroy(void);

double *readData(const char *filename, int N); // Reads the initial condition data [x y m vx vy L]
double *transform(const double *data, int N); // Allocates new 48N bytes in heap for the new order [m][x][y][vx][vy][L]
void SaveLastStep(const char *filename, double *DATA, int N); // Saves final [x y m vx vy L]
void SaveLastStep2(const char *filename, struct particle_s *current, double *L, int N);

void InitializeParticles(data_structure data_pointers);

#ifdef DEBUG
void printArray1D(double *v, int N);
void printAcceleration(int N);
void printPositions(int N);
#endif

int main(int argc, char *argv[]) {
    double start, end;
    long thread;
    pthread_t *thread_handles;
    N = atoi(argv[1]);
    const char *filename = argv[2];
    n_steps = atoi(argv[3]);
    dt = atof(argv[4]);
    G = 100.0 / N;
    M = atoi(argv[6]);

    double *data = readData(filename, N), *DATA = transform(data, N);
    double *m = DATA, *x = DATA + N, *y = DATA +2*N;
    double *u = DATA + 3*N, *v = DATA + 4*N;

    current = malloc(N * sizeof(struct particle_s));
    a = malloc(N * sizeof(vect_t));
    loc_a = malloc(M * N * sizeof(vect_t));
    thread_handles = malloc(M * sizeof(pthread_t));

    data_structure S = {m, x, y, u, v};
    InitializeParticles(S);

    Barrier_init();

    GET_TIME(start);
    for (thread = 0; thread < M; thread++)
        pthread_create(&thread_handles[thread], NULL, thread_function, (void *) thread);

    for (thread = 0; thread < M; thread++)
        pthread_join(thread_handles[thread], NULL);

    GET_TIME(end);

    Barrier_destroy();
    printf("Time elapsed: %lf\n", end - start);
    SaveLastStep2("../compare_gal_files/result_balanced.gal", current, DATA + 5*N, N);
    free(thread_handles);
    free(current);
    free(a);
    free(loc_a);
    free(data);
    free(DATA);
    return 0;
}

void Loop_schedule(int my_rank, int thread_count, int n, int sched, int *first_p, int *last_p, int *incr_p) {
    if (sched == CYCLIC) {
        *first_p = my_rank;
        *last_p = n;
        *incr_p = thread_count;
    } else {  /* sched == BLOCK */
        int quotient = n / thread_count;
        int remainder = n % thread_count;
        int my_iters;
        *incr_p = 1;
        if (my_rank < remainder) {
            my_iters = quotient + 1;
            *first_p = my_rank * my_iters;
        } else {
            my_iters = quotient;
            *first_p = my_rank * my_iters + remainder;
        }
        *last_p = *first_p + my_iters;
    }

}

void *thread_function(void *rank) {
    long my_rank = (long) rank;
    int thread;
    int n;
    int i;     /* Current particle               */
    int bfirst;   /* My first particle in blk sched */
    int blast;    /* My last particle in blk sched  */
    int bincr;    /* Loop increment in blk sched    */
    int cfirst;   /* My first particle in cyc sched */
    int clast;    /* My last particle in cyc sched  */
    int cincr;    /* Loop increment in cyc sched    */

    Loop_schedule(my_rank, M, N, BLOCK, &bfirst, &blast, &bincr);
    Loop_schedule(my_rank, M, N, CYCLIC, &cfirst, &clast, &cincr);

    for (n = 0; n < n_steps; n++) {
        memset(loc_a + my_rank * N, 0, N * sizeof(vect_t));
        Barrier();
        for (i = cfirst; i < clast; i += cincr)
            Compute_force(i, loc_a + my_rank * N);
        Barrier();
        for (i = bfirst; i < blast; i += bincr) {
            a[i][X] = a[i][Y] = 0.0;
            for (thread = 0; thread < M; thread++) {
                a[i][X] += loc_a[thread * N + i][X];
                a[i][Y] += loc_a[thread * N + i][Y];
            }
        }
        Barrier();
        for (i = bfirst; i < blast; i += bincr)
            Update_part(i);
        Barrier();
    }

    return NULL;
}





void Compute_force(int i, vect_t loc_a[]) {
    int j;
    vect_t f_ij;
    register double F1, F2;
    double R, R2, R3, dx2, dy2, invR3, Gx, Gy;
    double x_i = current[i].position[X];
    double y_i = current[i].position[Y];
    F1 = loc_a[i][X];
    F2 = loc_a[i][Y];
    double m_i = current[i].m;
#pragma GCC ivdep
    for (j = i + 1; j < N; j++) {
        f_ij[X] = current[j].position[X] - x_i;
        f_ij[Y] = current[j].position[Y] - y_i;
        //f_part_k[X] = dx; f_part_k[Y] = dy;
        dx2 = f_ij[X] * f_ij[X];
        dy2 = f_ij[Y] * f_ij[Y];
        R2 = dx2 + dy2;
        R = sqrt(R2) + eps;
        R3 = R * R;
        R3 *= R;
        invR3 = 1.0 / R3;
        Gx = Gy = invR3;
        f_ij[X] *= Gx;
        f_ij[Y] *= Gy;
        loc_a[j][X] -= f_ij[X] * m_i;
        loc_a[j][Y] -= f_ij[Y] * m_i;
        double m_j = current[j].m;
        F1 += f_ij[X] * m_j;
        F2 += f_ij[Y] * m_j;
    }
    loc_a[i][X] = F1;
    loc_a[i][Y] = F2;
}

void InitializeParticles(data_structure data_pointers) {
    int i;
    double *m = data_pointers.m;
    double *x = data_pointers.x;
    double *y = data_pointers.y;
    double *u = data_pointers.u;
    double *v = data_pointers.v;

    for (i = 0; i < N; i++) {
        current[i].m = m[i];
        current[i].position[X] = x[i];
        current[i].position[Y] = y[i];
        current[i].velocity[X] = u[i];
        current[i].velocity[Y] = v[i];
    }
}



void Update_part(int i) {
    current[i].velocity[X] += G * (a[i][X] * dt);
    current[i].position[X] += current[i].velocity[X] * dt;
    current[i].velocity[Y] += G * (a[i][Y] * dt);
    current[i].position[Y] += current[i].velocity[Y] * dt;
}

// BARRIER FUNCTIONALITY OK
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
void SaveLastStep2(const char *filename, struct particle_s *current, double *L, int N) {
    FILE *file = fopen(filename, "wb");
    if (file == NULL) {
        printf("Error opening file %s\n", filename);
        exit(1);
    }
    for (int i = 0; i < N; i++) {
        fwrite(&current[i].position[X], sizeof(double), 1, file);
        fwrite(&current[i].position[Y], sizeof(double), 1, file);
        fwrite(&current[i].m, sizeof(double), 1, file);
        fwrite(&current[i].velocity[X], sizeof(double), 1, file);
        fwrite(&current[i].velocity[Y], sizeof(double), 1, file);
        fwrite(&L[i], sizeof(double), 1, file);
    }
    fclose(file);
}

#ifdef DEBUG
void printArray1D(double *v, int N){
    for (int i = 0; i < N ; i++){
        printf("%lf\n", v[i]);
    }
}

void printAcceleration(int N){
    for (int i = 0; i < N; i ++)
        printf("%lf\n", a[i][X]);
}
void printPositions(int N){
    for (int i = 0; i < N; i++)
        printf("%lf\n", current[i].position[X]);
}
#endif