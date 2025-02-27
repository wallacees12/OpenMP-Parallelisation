//
// Created by Sam Wallace on 27.02.2025.
//
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

typedef struct{
  int start, end;
  int N;
  int n_threads;
  double *m;
  double *x;
  // Now each thread has a local acceleration array
  // double *a;
  double G;
  double eps;
  int thread;
  pthread_mutex_t *mutex;
} ThreadData;


static inline void
Compute_ax_ay(int i, int j,
              double *__restrict m,
              double *__restrict x, double *a, double G, double eps);

double *readData(const char *filename, int N);

double *transform(const double *data, int N);

void SaveLastStep(const char *filename, double *DATA, int N); // Doesn't step out of a mountain

void * compute_fores(void *arg){
  // Unroll all of the data
  ThreadData *data = (ThreadData *)arg;
  int start = data->start, end = data->end;
  int N = data->N;
  int n_threads = data->n_threads;
  double *m = data->m;
  double *x = data->x;
  // double *a = data->a; no longer needed as using local a
  double G = data->G;
  double eps = data->eps;
  int thread = data->thread;
  pthread_mutex_t *mutex = data->mutex;

  // Local acceleration
  double *local_a = calloc(2 * N, sizeof(double));

  // Thread will only work on its allocated particles in range (start,end)
  // printf("Thread started: start=%d, end=%d, N=%d\n", start, end, N);
  for(int i = 0; i < N-n_threads; i ++){
    int jchunk = (N - (i + 1) + n_threads - 1) / n_threads;
    int j_start = i + 1 + thread * jchunk;
    int j_end = i + 1 + (thread + 1) * jchunk;

    if (j_start >= N) continue;
    if (j_end > N) j_end = N;

    for(int j = j_start; j < j_end; j++){
        // Compute local i.e. ignores the mutex issues
        // Debugging prints
        // printf("Computing acceleration between i=%d, j=%d\n", i, j);
        // fflush(stdout); // Ensure prints appear before a crash
        Compute_ax_ay(i, j, m, x, local_a, G, eps);
    }
    // printf("Thread finished: start=%d, end=%d\n", start, end);
  }
  return local_a; // threaded function returns pointer, can be interpreted as VOID * we will cast
}



int main(int argc, char *argv[]) {

    if (argc != 7) {
        printf("Incorrect usage: ./gaslsim N filname nsteps delta_t graphics n_threads\n");
        return 1;}

    int N = atoi(argv[1]);
    const char *filename = argv[2];
    int n_steps = atoi(argv[3]);
    const double dt = atof(argv[4]), G = 100.0 / N, eps = 1e-3;
    // int graphics = atoi(argv[5]);
    int n_threads = atoi(argv[6]);

    // Sanity Check
    printf("n-threads = %d\n", n_threads);

    // Same as Vanilla
    double *data = readData(filename, N); // x0 y0 m0 vx0 vy0 L0
    double *DATA = transform(data, N); // [m0...mN-1] [x0 y0 x1 y1 ,....] [vx0 vy0, ...] [L0 L1,...]
    double *a = calloc(2 * N, sizeof(double)); // Now our Global acceleration
    double *m = DATA; // -> DATA[0] ; m[i]
    double *x = DATA + N; // -> x [x0, y0] -> x_i = x[2*i], y_i = x[2*i+1]
    double *v = DATA + 3 * N;

    // Initalise threads using pthreads
    pthread_t threads[n_threads];
    ThreadData thread_data[n_threads];
    // pthread_mutex_t mutex =  PTHREAD_MUTEX_INITIALIZER; No longer using this! We are calculating local then merge

    // Have to compute chunk size for each thread
    int chunk_size = N / n_threads;

    int i, j, n;
    n = 0;


    for(int n = 0; n < n_steps; n++) {

        for(int t=0; t < n_threads; t++){
            thread_data[t].start = t; // Start each thread at its index
            thread_data[t].end = N; // End at N particles
            thread_data[t].N = N;
            thread_data[t].m = m;
            thread_data[t].x = x;
            thread_data[t].G = G;
            thread_data[t].eps = eps;
            thread_data[t].n_threads = n_threads;
            thread_data[t].thread = t;
            pthread_create(&threads[t], NULL, compute_fores, &thread_data[t]);
        }

        memset(a, 0, 2 * N * sizeof(double));
        
        // Must wait for all our threads to come back to us for this step iteration
        
        for(int t = 0; t < n_threads; t++){
          double *local_a = NULL;
          pthread_join(threads[t], (void**)&local_a); // Retrieve our calculations
          for(int i = 0; i < 2 * N; i++){
            a[i] += local_a[i];
          }
          free(local_a);
        }
        double *accs = calloc(2*N,sizeof(double));
        for(int i=N-n_threads;i<N;i++){
            for(j=i+1; j<N;j++){
                Compute_ax_ay(i,j,m,x,accs,G,eps);
            }
            a[i]+=accs[i];  
        }
        free(accs);
        // We can definitely unroll this loop explicitly but the fucking compiler better
        // not be as lazy as me
        for (i = 0; i < 2 * N; i += 2) {
            v[i] += a[i] * dt; // Update x-velocity
            x[i] += v[i] * dt; // Update x-position
            v[i + 1] += a[i + 1] * dt; // Update y-velocity
            x[i + 1] += v[i + 1] * dt; // Update y-position
        }
        // Let's not even think about how this guy can be changed...
        // memset(a, 0, 2 * N * sizeof(double)); // At least we don't store accelerations in a matrix, that shit cringe
    }
    SaveLastStep("result_balanced_row.gal", DATA, N);
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
    dx = x[2 * j]     - x[2 * i];
    dy = x[2 * j + 1] - x[2 * i + 1];

    R2 = dx * dx + dy * dy;
    R = sqrt(R2) + eps;
    R3 = R * R * R;
    invR3 = 1.0 / R3;

    Gx = G * invR3 * dx;
    Gy = G * invR3 * dy;

    // For particle i, add contribution from j
    a[2 * i]     += Gx * m[j];
    a[2 * i + 1] += Gy * m[j];
    // For particle j, subtract contribution from i (Newton's third law)
    a[2 * j]     -= Gx * m[i];
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