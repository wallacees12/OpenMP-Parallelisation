#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

typedef struct{
  int start, end;
  int N;
  double *m;
  double *x;
  double *a;
  double G;
  double eps;
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
  double *m = data->m;
  double *x = data->x;
  double *a = data->a;
  double G = data->G;
  double eps = data->eps;
  pthread_mutex_t *mutex = data->mutex;

  // Thread will only work on its allocated particles in range (start,end)
  printf("Thread started: start=%d, end=%d, N=%d\n", start, end, N);
  for(int i = start; i < end; i++){
    for(int j = i + 1; j < N; j++){
        // Compute local i.e. ignores the mutex issues
        // Debugging prints
        printf("Computing acceleration between i=%d, j=%d\n", i, j);
        fflush(stdout); // Ensure prints appear before a crash

        double local_a[2] = {0.0, 0.0};
        Compute_ax_ay(i, j, m, x, local_a, G, eps);

        // Ensure indices are valid before accessing a[]
        if (2 * i >= 2 * N || 2 * j >= 2 * N) {
            printf("Error: Out-of-bounds access! i=%d, j=%d\n", i, j);
            exit(1);
        }
        pthread_mutex_lock(mutex);
        a[2 * i] += local_a[0];
        a[2 * i + 1] += local_a[1];
        a[2 * j] -= local_a[0];
        a[2 * j + 1] -= local_a[1];
        pthread_mutex_unlock(mutex);
    }
    printf("Thread finished: start=%d, end=%d\n", start, end);
  }
  return NULL; // threaded function must return NULL
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
    double *a = calloc(2 * N, sizeof(double)); // [ax0 ay0, ....]
    double *m = DATA; // -> DATA[0] ; m[i]
    double *x = DATA + N; // -> x [x0, y0] -> x_i = x[2*i], y_i = x[2*i+1]
    double *v = DATA + 3 * N;

    // Initalise threads using pthreads
    pthread_t threads[n_threads];
    ThreadData thread_data[n_threads];
    pthread_mutex_t mutex =  PTHREAD_MUTEX_INITIALIZER;

    // Have to compute chunk size for each thread
    int chunk_size = N / n_threads;

    int i, j, n;
    n = 0;
    // Change from a while loop to a for loop @juan
    for(int n = 0; n < n_steps; n++) {
      printf("Iteration %d\n", n);
        // Threading time
        for(int t=0l; t < n_threads; t++){
          thread_data[t].start = t * chunk_size;
          // Fancy way to write if it is last threads go to end N or else calculate end
          thread_data[t].end = (t == n_threads - 1) ? N : (t + 1) * chunk_size;
          thread_data[t].N = N;
          thread_data[t].m = m;
          thread_data[t].x = x;
          thread_data[t].a = a;
          thread_data[t].G = G;
          thread_data[t].eps = eps;
          thread_data[t].mutex = &mutex;
          pthread_create(&threads[t], NULL, compute_fores, &thread_data[t]);
          //BOOM THREAD MADE
          printf("Thread %d is live\n", t+1);
        }

        // Must wait for all our threads to come back to us for this step iteration
        for(int t = 0; t < n_threads; t++){
          pthread_join(threads[t], NULL);
        }

        // We can definitely unroll this loop explicitly but the fucking compiler better
        // not be as lazy as me
        for (i = 0; i < 2 * N; i += 2) {
            v[i] += a[i] * dt; // Update x-velocity
            x[i] += v[i] * dt; // Update x-position
            v[i + 1] += a[i + 1] * dt; // Update y-velocity
            x[i + 1] += v[i + 1] * dt; // Update y-position
        }
        // Let's not even think about how this guy can be changed...
        memset(a, 0, 2 * N * sizeof(double)); // At least we don't store accelerations in a matrix, that shit cringe
    }
    SaveLastStep("result.gal", DATA, N);
    free(a);
    free(data);
    free(DATA);

    // Only update is detroy the mutex
    pthread_mutex_destroy(&mutex);
    return 0;
}

static inline void
Compute_ax_ay(int i, int j, double *__restrict m, double *__restrict x, double *a, double G, double eps) {
    int I, J;


    double dx, dx2, dy, dy2, invR3, Gx, Gy, R, R2, R3, m_i = m[i], m_j = m[j];
    I = 2 * i;
    J = 2 * j;

    // Debugging: Check if indices are valid
    if (I >= 2 * 1000 || J >= 2 * 1000) {  // Assuming N=1000 for debug
        printf("ERROR: Out-of-bounds access in Compute_ax_ay: i=%d, j=%d, I=%d, J=%d\n", i, j, I, J);
        exit(1);
    }

    dx = x[J] - x[I];
    dx2 = dx * dx;
    dy = x[J + 1] - x[I + 1];
    dy2 = dy * dy;
    R2 = dx2 + dy2;
    R = sqrt(R2) + eps;
    R3 = R * R;
    R3 *= R;
    invR3 = 1.0 / R3;
    Gx = Gy = G * invR3; // Theoretically G should be able to be put in the update loop, but it no work
    Gx *= dx;
    Gy *= dy;

    a[I] += Gx * m_j;
    a[J] -= Gx * m_i;
    a[I + 1] += Gy * m_j;
    a[J + 1] -= Gy * m_i;


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