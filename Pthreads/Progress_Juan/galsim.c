//
// Created by Juan Rodriguez on 2025-02-18.
//

/* GRAPHICS OMITTED */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "timer.h"

double *readData(const char *filename, int N); // Reads the initial condition data [x y m vx vy L]
double *transform(const double *data, int N); // Allocates new 48N bytes in heap for the new order [m][x][y][vx][vy][L]
void SaveLastStep(const char *filename, double *DATA, int N); // Saves final [x y m vx vy L]

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
    int graphics = atoi(argv[5]);

    double *data = readData(filename, N);
    double *DATA = transform(data, N);
    // [m][x][y][u][v] ~ [ax][ay]
    double *a = calloc(2 * N, sizeof(double)); // Allocate a different buffer for accelerations i x and y
    double *m = DATA; // Pointer to first mass in DATA = [m][x][y][vx][vy][L]
    double *x = DATA + N; // Pointer to first x coordinate  in DATA = [m][x][y][vx][vy][L]
    double *y = DATA + 2 * N; // Pointer to first y coordinate  in DATA = [m][x][y][vx][vy][L]
    double *u = DATA + 3 * N;  // Pointer to first x velocity  in DATA = [m][x][y][vx][vy][L]
    double *v = DATA + 4 * N; // Pointer to first y velocity  in DATA = [m][x][y][vx][vy][L]
    double *ax = a; // Pointer to first x acceleration in a = [ax][ay]
    double *ay = a + N;// Pointer to first y acceleration in a = [ax][ay]

    double m_j, x_i, y_i;
    int i, j, n;
    // Define frequently used variables as register variables in case of the compiler doesn't do this automatically
    register double F1, F2, dx, dy, dx2, dy2, R2, R, R3, invR3, Gx, Gy, m_i; // Frequently used in innermost loop
    n = 0;
    double start, end;
    GET_TIME(start);
    while (n < n_steps) {
        // Index tuples will be ordered as an upper triangular matrix representation
        for (i = 0; i < N; i++) {
            x_i = x[i]; y_i = y[i]; m_i = m[i]; // Auto-vectorization works better with i,j independence
            F1 = ax[i]; F2 = ay[i];
#pragma GCC ivdep
            for (j = i + 1; j < N; j++) { // Iterations in j decreases as N-1, N-2, N-3, ..., 3,  2, 1
                // Compute quantities need for accelerations in x and y
                dx = x[j] - x_i; dy = y[j] - y_i;
                dx2 = dx * dx; dy2 = dy * dy;
                R2 = dx2 + dy2;
                R = sqrt(R2) + eps;
                R3 = R * R;
                R3 *= R;
                invR3 = 1.0 / R3;
                Gx = Gy = invR3;
                Gx *= dx; Gy *= dy;
                // Sum over the accelerations on particle i from particle j
                ax[j] -= Gx * m_i; ay[j] -= Gy * m_i;
                m_j = m[j];
                // Compute the equal in magnitude but opposite directed forces (skew symmetry)
                F1 += Gx * m_j; F2 += Gy * m_j; // F1 and F2 are accumulators that make the loops independent
            }
            // Symplectic Euler
            u[i] += G * (F1 * dt); // Update x velocity, F1, F2, should be read as accelerations not forces
            x[i] += u[i] * dt; // Update x position
            v[i] += G * (F2 * dt); // Update y velocity
            y[i] += v[i] * dt; // Update y position
        }
        n++; // Increment time step counter
        memset(a, 0, 2 * N * sizeof(double)); // Set acceleration buffer elements to 0
    }
    GET_TIME(end);
    printf("Time ellapsed: %lf\n", end-start);
    SaveLastStep("result.gal", DATA, N);
    free(a);
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