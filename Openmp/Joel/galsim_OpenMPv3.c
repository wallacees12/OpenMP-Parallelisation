#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>

//Pretty much Sams code

typedef struct{
    int N;
    double *m;
    double *x;
    double *y;
    double G;
    double eps;
} ThreadedData;

static double get_wall_seconds() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
    return seconds;
}

static inline void Compute_ax_ay(int i, int j, double m_i, double m_j,
                                 double  x_i, double  x_j, 
                                 double y_i, double y_j,
                                 double *ax_i, double *ay_i,
                                 double *ax_j, double *ay_j, double G, double eps);

double *readData(const char *filename, int N);
double *transform(const double *data, int N);
void SaveLastStep(const char *filename, double *DATA, int N);

void compute_forces(ThreadedData *data, double *a, int n_threads) {
    int N = data->N;
    double *m = data->m;
    double *x = data->x;
    double *y = data->y;
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
            double m_i=m[i]; double x_i=x[i]; double y_i=y[i];
            for (int j = i + 1; j < N; j++) {
                double m_j=m[j]; double x_j=x[j]; double y_j=y[j];
                double ax_i = 0.0, ay_i = 0.0, ax_j = 0.0, ay_j = 0.0;
                Compute_ax_ay(i, j, m_i, m_j, x_i, x_j, y_i, y_j, &ax_i, &ay_i, &ax_j, &ay_j, G, eps);

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
    double *y = DATA + 2 * N;
    double *u = DATA + 3 * N;
    double *v = DATA + 4 * N;

    ThreadedData threadedData = {N, m, x, y, G, eps};
    double start = get_wall_seconds();

    for (int n = 0; n < n_steps; n++) {
        memset(a, 0, 2 * N * sizeof(double));
        compute_forces(&threadedData, a, n_threads);

        #pragma omp parallel for num_threads(n_threads)
        for (int i = 0; i < N; i ++) {
            v[i] += a[i] * dt;
            x[i] += v[i] * dt;
            u[i] += a[i + 1] * dt;
            y[i] += v[i] * dt;
        }
    }
    printf("Time taken: %.3lfs\n", get_wall_seconds()-start);
    SaveLastStep("galsim_OpenMP.gal", DATA, N);
    free(a);
    free(data);
    free(DATA);
    return 0;
}

static inline void Compute_ax_ay(int i, int j, double  m_i, double  m_j,
                                 double  x_i, double  x_j, 
                                 double  y_i, double  y_j,
                                 double *ax_i, double *ay_i,
                                 double *ax_j, double *ay_j, double G, double eps){
    double dx = x_j - x_i;
    double dy = y_j - y_i;

    double R2 = dx * dx + dy * dy;
    double R = sqrt(R2) + eps;
    double R3 = R * R * R;
    double invR3 = 1.0 / R3;

    double Gx = G * invR3 * dx;
    double Gy = G * invR3 * dy;

    *ax_i = Gx * m_j;
    *ay_i = Gy * m_j;
    *ax_j = -Gx * m_i;
    *ay_j = -Gy * m_i;
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