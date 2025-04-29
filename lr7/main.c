#include <stdio.h>
#include <math.h>
#include <gsl/gsl_interp.h>

#define N 11

// Исходные данные
double x[N] = {0.7416, 0.9752, 1.1664, 1.5400, 3.7888, 3.8776, 4.5128, 4.9328, 5.1840, 5.5352, 5.6976};
double y[N] = {-0.7344, -0.0515, 0.2543, 0.3627, 3.9350, 4.7417, 13.3534, 22.2931, 29.0949, 40.6465, 46.8569};
double target_x = 5.4008;

// Интерполяция Лагранжа
double lagrange_interpolation(double xi) {
    double result = 0.0;
    for (int i = 0; i < N; i++) {
        double term = y[i];
        for (int j = 0; j < N; j++) {
            if (j != i) {
                term *= (xi - x[j]) / (x[i] - x[j]);
            }
        }
        result += term;
    }
    return result;
}

// Интерполяция Ньютона (разделенные разности)
double newton_interpolation(double xi) {
    double diff_table[N][N];
    
    for (int i = 0; i < N; i++) {
        diff_table[i][0] = y[i];
    }
    
    for (int j = 1; j < N; j++) {
        for (int i = 0; i < N - j; i++) {
            diff_table[i][j] = (diff_table[i+1][j-1] - diff_table[i][j-1]) / (x[i+j] - x[i]);
        }
    }
    
    double result = diff_table[0][0];
    double product = 1.0;
    
    for (int i = 1; i < N; i++) {
        product *= (xi - x[i-1]);
        result += diff_table[0][i] * product;
    }
    
    return result;
}

// Интерполяция Эйткена
double aitken_interpolation(double xi) {
    double P[N][N];
    
    for (int i = 0; i < N; i++) {
        P[i][0] = y[i];
    }
    
    for (int k = 1; k < N; k++) {
        for (int i = 0; i < N - k; i++) {
            P[i][k] = ((xi - x[i]) * P[i+1][k-1] - (xi - x[i+k]) * P[i][k-1]) / (x[i+k] - x[i]);
        }
    }
    
    return P[0][N-1];
}

// Полиномиальная интерполяция через GSL
double gsl_polynomial_interpolation(double xi) {
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_interp *interp = gsl_interp_alloc(gsl_interp_polynomial, N);
    
    gsl_interp_init(interp, x, y, N);
    double result = gsl_interp_eval(interp, x, y, xi, acc);
    
    gsl_interp_free(interp);
    gsl_interp_accel_free(acc);
    
    return result;
}

void print_comparison(const char* method_name, double value, double gsl_value) {
    double diff = fabs(value - gsl_value);
    printf("%-15s: %.16f  (diff with GSL: %f)\n", method_name, value, diff);
}

int main() {
	for (int i = 40; i < 50; i++) {
		double target_x = i * 0.1;
    		double gsl_result = gsl_polynomial_interpolation(target_x);
    		printf("Interpolation at x = %.4f\n", target_x);
    		print_comparison("Lagrange", lagrange_interpolation(target_x), gsl_result);
    		print_comparison("Newton", newton_interpolation(target_x), gsl_result);
   	 	print_comparison("Aitken", aitken_interpolation(target_x), gsl_result);
    		print_comparison("GSL", gsl_result, gsl_result);
	}
    return 0;
}
