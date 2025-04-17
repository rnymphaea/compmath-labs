// gaussian_jacobi_solver.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#define MAX_SIZE 300

int solve_gauss(float matrix[MAX_SIZE][MAX_SIZE], float vector[MAX_SIZE], float solution[MAX_SIZE], int size) {
    int i, j, k, max_row = 0;
    int error = 0;
    float max_val, temp, factor, sum;

    if (size < 2) { error = 2; return error; }
    if (size > MAX_SIZE) { error = 3; return error; }

    for (k = 0; k < size - 1; k++) {
        max_val = 0.0;
        for (i = k; i < size; i++) {
            if (fabs(matrix[i][k]) > fabs(max_val)) {
                max_val = matrix[i][k];
                max_row = i;
            }
        }
        if (max_val == 0.0) { error = 1; return error; }

        if (max_row != k) {
            for (j = k; j < size; j++) {
                temp = matrix[max_row][j];
                matrix[max_row][j] = matrix[k][j];
                matrix[k][j] = temp;
            }
            temp = vector[max_row];
            vector[max_row] = vector[k];
            vector[k] = temp;
        }

        for (i = k + 1; i < size; i++) {
            factor = -matrix[i][k] / matrix[k][k];
            matrix[i][k] = 0.0;
            for (j = k + 1; j < size; j++) {
                matrix[i][j] += matrix[k][j] * factor;
            }
            vector[i] += vector[k] * factor;
        }
    }

    if (matrix[size - 1][size - 1] == 0.0) { error = 1; return error; }

    for (i = size - 1; i >= 0; i--) {
        sum = 0.0;
        for (j = i + 1; j < size; j++) {
            sum += matrix[i][j] * solution[j];
        }
        solution[i] = (vector[i] - sum) / matrix[i][i];
    }

    return error;
}

int solve_iterative(float A[MAX_SIZE][MAX_SIZE], float b[MAX_SIZE], float x[MAX_SIZE], int n, float epsilon, int *iterations) {
    int i, j, error = 0;
    float prev_x[MAX_SIZE];
    float diff;
    double norm;

    if (n < 2) { error = 2; return error; }
    if (n > MAX_SIZE) { error = 3; return error; }

    gsl_matrix *temp = gsl_matrix_alloc(n, n);
    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++){
            gsl_matrix_set(temp, i, j, (double)A[i][j]);
        }
    }
    norm = gsl_matrix_norm1(temp);
    gsl_matrix_free(temp);

    if (norm >= 1.0) { error = 1; return error; }

    for (i = 0; i < n; i++) {
        x[i] = b[i];
    }
    *iterations = 0;

    do {
        for (i = 0; i < n; i++) {
            prev_x[i] = x[i];
        }

        for (i = 0; i < n; i++) {
            float sum = 0.0f;
            for (j = 0; j < n; j++) {
                if (j != i)
                    sum += A[i][j] * prev_x[j];
            }
            x[i] = (b[i] - sum) / A[i][i];
        }

        (*iterations)++;

        gsl_vector *vec_diff = gsl_vector_alloc(n);
        for (i = 0; i < n; i++) {
            gsl_vector_set(vec_diff, i, (double)(x[i] - prev_x[i]));
        }
        diff = gsl_blas_dnrm2(vec_diff);
        gsl_vector_free(vec_diff);

    } while (diff > epsilon);

    return error;
}

int main(void) {
    int size, i, j, result, choice;
    float matrix[MAX_SIZE][MAX_SIZE], vector[MAX_SIZE];
    float iter_matrix[MAX_SIZE][MAX_SIZE], iter_vector[MAX_SIZE];
    float gauss_solution[MAX_SIZE], iter_solution[MAX_SIZE];
    double gsl_matrix_d[MAX_SIZE][MAX_SIZE], gsl_vector_d[MAX_SIZE];

    gsl_matrix *A_gsl = NULL, *A_original = NULL, *A_inverse = NULL;
    gsl_vector *b_gsl = NULL, *x_gsl = NULL, *difference = NULL, *diff_vec = NULL;
    gsl_permutation *perm = NULL;

    const gsl_rng_type *rng_type;
    gsl_rng *rng;
    gsl_rng_env_setup();
    rng_type = gsl_rng_default;
    rng = gsl_rng_alloc(rng_type);
    gsl_rng_set(rng, time(NULL));

    printf("Enter system size (n): ");
    if (scanf("%d", &size) != 1 || size < 2 || size > MAX_SIZE) {
        printf("Invalid size.\n");
        gsl_rng_free(rng);
        return -1;
    }

    printf("Select input method:\n1. Manual input\n2. Random generation\n");
    if (scanf("%d", &choice) != 1) {
        printf("Input error.\n");
        return -1;
    }

    if (choice == 1) {
        for(i = 0; i < size; i++)
            for(j = 0; j < size; j++) {
                scanf("%f", &matrix[i][j]);
                iter_matrix[i][j] = matrix[i][j];
                gsl_matrix_d[i][j] = matrix[i][j];
            }
        for(i = 0; i < size; i++) {
            scanf("%f", &vector[i]);
            iter_vector[i] = vector[i];
            gsl_vector_d[i] = vector[i];
        }
    } else if (choice == 2) {
        for (i = 0; i < size; i++) {
            float row_sum = 0.0f;
            for (j = 0; j < size; j++) {
                if (i != j) {
                    float val = (gsl_rng_uniform(rng) * 10.0f - 5.0f) / 10000;
                    matrix[i][j] = val;
                    iter_matrix[i][j] = val;
                    gsl_matrix_d[i][j] = val;
                    row_sum += fabs(val);
                }
            }
            // Диагональное преобладание
            float diag = (row_sum + gsl_rng_uniform(rng) * 5.0f + 1.0f) / 1000;
            matrix[i][i] = diag;
            iter_matrix[i][i] = diag;
            gsl_matrix_d[i][i] = diag;

            vector[i] = gsl_rng_uniform(rng) * 20.0f - 10.0f;
            iter_vector[i] = vector[i];
            gsl_vector_d[i] = vector[i];
        }

        printf("Generated matrix A and vector b:\n");
        for(i = 0; i < size; i++) {
            for(j = 0; j < size; j++)
                printf("%8.3f ", matrix[i][j]);
            printf("| %8.3f\n", vector[i]);
        }
    } else {
        printf("Invalid choice.\n");
        return -1;
    }

    A_gsl = gsl_matrix_alloc(size, size);
    b_gsl = gsl_vector_alloc(size);
    x_gsl = gsl_vector_alloc(size);
    perm = gsl_permutation_alloc(size);
    for(i = 0; i < size; i++) {
        for(j = 0; j < size; j++)
            gsl_matrix_set(A_gsl, i, j, gsl_matrix_d[i][j]);
        gsl_vector_set(b_gsl, i, gsl_vector_d[i]);
    }

    int signum;
    gsl_linalg_LU_decomp(A_gsl, perm, &signum);
    gsl_linalg_LU_solve(A_gsl, perm, b_gsl, x_gsl);

    printf("\nGSL solution:\n");
    for(i = 0; i < size; i++)
        printf("x_gsl[%d] = %.20f\n", i, gsl_vector_get(x_gsl, i));

    A_original = gsl_matrix_alloc(size, size);
    for(i = 0; i < size; i++)
        for(j = 0; j < size; j++)
            gsl_matrix_set(A_original, i, j, gsl_matrix_d[i][j]);

    double matrix_norm = gsl_matrix_norm1(A_original);
    A_inverse = gsl_matrix_alloc(size, size);
    if (gsl_linalg_LU_invert(A_gsl, perm, A_inverse) != 0)
        printf("Failed to compute inverse matrix.\n");

    double inverse_norm = gsl_matrix_norm1(A_inverse);
    printf("\nMatrix condition number: %f\n", matrix_norm * inverse_norm);

    result = solve_gauss(matrix, vector, gauss_solution, size);
    if (result == 0) {
        printf("\nGauss method solution:\n");
        for(i = 0; i < size; i++)
            printf("x_gauss[%d] = %.20f\n", i, gauss_solution[i]);
    } else {
        printf("Error solving using Gauss method (code: %d)\n", result);
        return -1;
    }

    difference = gsl_vector_alloc(size);
    printf("\nDifference between Gauss and GSL:\n");
    for(i = 0; i < size; i++) {
        double diff = gauss_solution[i] - gsl_vector_get(x_gsl, i);
        printf("x[%d]: %.20f\n", i, diff);
        gsl_vector_set(difference, i, diff);
    }
    printf("Norm: %.20f\n", gsl_blas_dnrm2(difference));
    gsl_vector_free(difference);

    float epsilon;
    int iterations;
    printf("\nEnter required precision (epsilon) for iterative method: ");
    scanf("%f", &epsilon);
    result = solve_iterative(iter_matrix, iter_vector, iter_solution, size, epsilon, &iterations);
    if(result == 0) {
        printf("\nIterative method solution:\n");
        for(i = 0; i < size; i++)
            printf("x_iter[%d] = %.20f\n", i, iter_solution[i]);
        printf("Iterations: %d\n", iterations);
    } else {
        printf("Iterative method failed (code: %d)\n", result);
    }

    float epsilon_values[] = {0.1f, 0.01f, 0.001f, 0.0001f, 0.00001f, 0.000001f};
    printf("\nEpsilon vs Iterations and Error:\n");
    for(i = 0; i < 6; i++) {
        float temp_solution[MAX_SIZE];
        int iters;
        result = solve_iterative(iter_matrix, iter_vector, temp_solution, size, epsilon_values[i], &iters);
        if (result != 0) continue;
        diff_vec = gsl_vector_alloc(size);
        for(j = 0; j < size; j++) {
            gsl_vector_set(diff_vec, j, temp_solution[j] - gsl_vector_get(x_gsl, j));
        }
        double diff_norm = gsl_blas_dnrm2(diff_vec);
        gsl_vector_free(diff_vec);
        printf("epsilon = %.6f, iterations = %d, error norm = %.20f\n", epsilon_values[i], iters, diff_norm);
    }

    gsl_matrix_free(A_gsl);
    gsl_vector_free(b_gsl);
    gsl_vector_free(x_gsl);
    gsl_permutation_free(perm);
    gsl_matrix_free(A_original);
    gsl_matrix_free(A_inverse);
    gsl_rng_free(rng);

    return 0;
}

