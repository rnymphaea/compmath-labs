#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#define MAX_SIZE 700  // Максимальный размер системы
#define EPS 0.0000001

// Реализация метода Гаусса с выбором главного элемента
int solve_gaussian_elimination(double matrix[MAX_SIZE][MAX_SIZE], 
                              double vector[MAX_SIZE], 
                              double solution[MAX_SIZE], 
                              int system_size) {
    int row, col, step, max_row = 0;
    int error_flag = 0;  // Флаг ошибки (0 - успех)

    // Проверка допустимости размера системы
    if (system_size < 2) return 2; 
    if (system_size > MAX_SIZE) return 3;

    // Прямой ход метода Гаусса
    for (step = 0; step < system_size - 1; step++) {
        // Поиск строки с максимальным элементом
        double max_val = 0.0;
        for (row = step; row < system_size; row++) {
            if (fabs(matrix[row][step]) > fabs(max_val)) {
                max_val = matrix[row][step];
                max_row = row;
            }
        }
        
        // Проверка на вырожденность
        if (fabs(max_val - 0) < EPS) return 1;
        
        // Перестановка строк при необходимости
        if (max_row != step) {
            for (col = step; col < system_size; col++) {
                double temp = matrix[max_row][col];
                matrix[max_row][col] = matrix[step][col];
                matrix[step][col] = temp;
            }
            double temp = vector[max_row];
            vector[max_row] = vector[step];
            vector[step] = temp;
        }

        // Исключение переменной
        for (row = step + 1; row < system_size; row++) {
            double factor = -matrix[row][step] / matrix[step][step];
            matrix[row][step] = 0.0;
            for (col = step + 1; col < system_size; col++) {
                matrix[row][col] += matrix[step][col] * factor;
            }
            vector[row] += vector[step] * factor;
        }
    }

    // Проверка последнего диагонального элемента
    if (matrix[system_size-1][system_size-1] == 0.0) return 1;

    // Обратный ход
    for (row = system_size - 1; row >= 0; row--) {
        double sum = 0.0;
        for (col = row + 1; col < system_size; col++) {
            sum += matrix[row][col] * solution[col];
        }
        solution[row] = (vector[row] - sum) / matrix[row][row];
    }

    return error_flag;
}

// Вычисление нормы матрицы (максимальная сумма по строкам)
double compute_matrix_norm(const gsl_matrix *mat) {
    double max_norm = 0.0;
    for (size_t i = 0; i < mat->size1; i++) {
        double row_sum = 0.0;
        for (size_t j = 0; j < mat->size2; j++) {
            row_sum += fabs(gsl_matrix_get(mat, i, j));
        }
        max_norm = row_sum > max_norm ? row_sum : max_norm;
    }
    return max_norm;
}

// Функция для вывода матрицы и вектора
void print_system(double matrix[MAX_SIZE][MAX_SIZE], double vector[MAX_SIZE], int size) {
    printf("\nВведенная система:\n");
    for (int i = 0; i < size; i++) {
        printf("| ");
        for (int j = 0; j < size; j++) {
            printf("%8.3f ", matrix[i][j]);
        }

        printf("| | x%-2d | = | %8.3f |\n", i+1, vector[i]);
    }
}

int main() {
    int system_size, i, j, option;
    double matrix[MAX_SIZE][MAX_SIZE], vector[MAX_SIZE], solution[MAX_SIZE];

    // Инициализация генератора случайных чисел
    const gsl_rng_type *rng_type;
    gsl_rng *rng;
    gsl_rng_env_setup();
    rng_type = gsl_rng_default;
    rng = gsl_rng_alloc(rng_type);
    
    gsl_rng_set(rng, time(NULL));
    printf("Введите размер системы (2-%d): ", MAX_SIZE);
    if (scanf("%d", &system_size) != 1 || system_size < 2 || system_size > MAX_SIZE) {
        printf("Некорректный размер системы\n");
        gsl_rng_free(rng);
        return EXIT_FAILURE;
    }

    printf("Выберите способ ввода:\n1. Ручной\n2. Случайный\n> ");
    if (scanf("%d", &option) != 1 || (option != 1 && option != 2)) {
        printf("Некорректный выбор\n");
        gsl_rng_free(rng);
        return EXIT_FAILURE;
    }

    // Ввод/генерация матрицы и вектора
    if (option == 1) {
        printf("Введите матрицу %dx%d:\n", system_size, system_size);
        for (i = 0; i < system_size; i++) {
            for (j = 0; j < system_size; j++) {
                scanf("%lf", &matrix[i][j]);
            }
        }
        printf("Введите вектор правой части:\n");
        for (i = 0; i < system_size; i++) {
            scanf("%lf", &vector[i]);
        }
    } else {
        for (i = 0; i < system_size; i++) {
            for (j = 0; j < system_size; j++) {
                matrix[i][j] = gsl_rng_uniform(rng) * 20.0 - 10.0;
            }
            vector[i] = gsl_rng_uniform(rng) * 20.0 - 10.0;
        }
    }

    // Вывод введенной системы
    print_system(matrix, vector, system_size);

    // Решение системы методом Гаусса
    int status = solve_gaussian_elimination(matrix, vector, solution, system_size);
    if (status != 0) {
        printf("Ошибка %d при решении системы\n", status);
        gsl_rng_free(rng);
        return EXIT_FAILURE;
    }

    printf("\nРешение системы (метод Гаусса):\n");
    for (i = 0; i < system_size; i++) {
        printf("x[%d] = %f\n", i, solution[i]);
    }

    // Решение через GSL
    gsl_matrix *gsl_mat = gsl_matrix_alloc(system_size, system_size);
    gsl_vector *gsl_vec = gsl_vector_alloc(system_size);
    gsl_vector *gsl_sol = gsl_vector_alloc(system_size);
    gsl_permutation *perm = gsl_permutation_alloc(system_size);
    int signum;

    for (i = 0; i < system_size; i++) {
        for (j = 0; j < system_size; j++) {
            gsl_matrix_set(gsl_mat, i, j, matrix[i][j]);
        }
        gsl_vector_set(gsl_vec, i, vector[i]);
    }

    gsl_linalg_LU_decomp(gsl_mat, perm, &signum);
    gsl_linalg_LU_solve(gsl_mat, perm, gsl_vec, gsl_sol);

    printf("\nРешение (GSL):\n");
    for (i = 0; i < system_size; i++) {
        printf("x[%d] = %g\n", i, gsl_vector_get(gsl_sol, i));
    }
    
     // Вычисление разницы между решениями и нормы разницы (евклидова норма)
    double diff_norm = 0.0;
    printf("\nРазница решений (метод Гаусса - GSL):\n");
    for(i = 0; i < system_size; i++){
        double diff = solution[i] - gsl_vector_get(gsl_sol, i);
        printf("Разница для x[%d] = %.20f\n", i, diff);
        diff_norm += diff * diff;
    }
    diff_norm = sqrt(diff_norm);
    printf("\nНорма разницы: %.20f\n", diff_norm);

    // Вычисление числа обусловленности
    gsl_matrix *original_mat = gsl_matrix_alloc(system_size, system_size);
    for (i = 0; i < system_size; i++) {
        for (j = 0; j < system_size; j++) {
            gsl_matrix_set(original_mat, i, j, matrix[i][j]);
        }
    }

    gsl_matrix *inverse_mat = gsl_matrix_alloc(system_size, system_size);
    gsl_linalg_LU_invert(gsl_mat, perm, inverse_mat);

    double mat_norm = compute_matrix_norm(original_mat);
    double inv_norm = compute_matrix_norm(inverse_mat);
    double cond_number = mat_norm * inv_norm;

    printf("\nЧисло обусловленности: %f\n", cond_number);

    // Освобождение ресурсов
    gsl_matrix_free(gsl_mat);
    gsl_matrix_free(original_mat);
    gsl_matrix_free(inverse_mat);
    gsl_vector_free(gsl_vec);
    gsl_vector_free(gsl_sol);
    gsl_permutation_free(perm);
    gsl_rng_free(rng);

    return EXIT_SUCCESS;
}
