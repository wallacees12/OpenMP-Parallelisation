//
// Created by Juan Rodriguez on 2025-02-14.
//

//
// Created by Juan Rodriguez on 2025-02-13.
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


static inline void
Compute_ax_ay(int i, int j,
              double *__restrict m,
              double *__restrict x, double *a, double G, double eps);

double *readData(const char *filename, int N);

double *transform(const double *data, int N);

void SaveLastStep(const char *filename, double *DATA, int N); // Doesn't step out of a mountain


int main(int argc, char *argv[]) {

    if (argc < 6) {
        printf("Too few of arguments.\n");
        return 1;
    } else if (6 < argc) {
        printf("Too many arguments.\n");
        return 1;
    }

    int N = atoi(argv[1]);
    const char *filename = argv[2];
    int n_steps = atoi(argv[3]);
    const double dt = atof(argv[4]), G = 100.0 / N, eps = 1e-3;
    // int graphics = atoi(argv[5]);

    double *data = readData(filename, N); // x0 y0 m0 vx0 vy0 L0
    double *DATA = transform(data, N); // [m0...mN-1] [x0 y0 x1 y1 ,....] [vx0 vy0, ...] [L0 L1,...]
    double *a = calloc(2 * N, sizeof(double)); // [ax0 ay0, ....]
    double *m = DATA; // -> DATA[0] ; m[i]
    double *x = DATA + N; // -> x [x0, y0] -> x_i = x[2*i], y_i = x[2*i+1]
    double *v = DATA + 3 * N;

    int i, j, n;
    n = 0;
    while (n < n_steps) {
        for (i = 0; i < N; i++) {
            // Might need to do a jailbreak for x_i, y_i, m_i? ANSWER: NO compiler handles it.
            for (j = i + 1; j < N; j++) {
                // This function shit is also inlined, it has a fat dad body. 3 double pointers too many?
                // Well fuck m and x are restricted for hopefully achieve some pointer aliasing
                // since after all m and x are stored contiguously. [m][xy|xy|...|xy]
                Compute_ax_ay(i, j, m, x, a, G, eps); // Inlined function as efficient as not inline function
            }
        }
        // We can definitely unroll this loop explicitly but the fucking compiler better
        // not be as lazy as me
        for (i = 0; i < 2 * N; i += 2) {
            v[i] += a[i] * dt; // Update x-velocity
            x[i] += v[i] * dt; // Update x-position
            v[i + 1] += a[i + 1] * dt; // Update y-velocity
            x[i + 1] += v[i + 1] * dt; // Update y-position
        }
        n++;
        // Let's not even think about how this guy can be changed...
        memset(a, 0, 2 * N * sizeof(double)); // At least we don't store accelerations in a matrix, that shit cringe
    }
    SaveLastStep("result.gal", DATA, N);
    free(a);
    free(data);
    free(DATA);
    return 0;
}

static inline void
Compute_ax_ay(int i, int j, double *__restrict m, double *__restrict x, double *a, double G, double eps) {
    int I, J;
    double dx, dx2, dy, dy2, invR3, Gx, Gy, R, R2, R3, m_i = m[i], m_j = m[j];
    I = 2 * i;
    J = 2 * j;
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