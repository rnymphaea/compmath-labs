#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double f(double x) {
    return cos(x * x) - 10.0 * x;
}

void computeDividedDifferences(int n, double x[], double y[], double coeffs[]) {
    for (int i = 0; i < n; i++) {
        coeffs[i] = y[i];
    }
    for (int j = 1; j < n; j++) {
        for (int i = n - 1; i >= j; i--) {
            coeffs[i] = (coeffs[i] - coeffs[i - 1]) / (x[i] - x[i - j]);
        }
    }
}

double newtonInterpolate(double x, int n, double nodes[], double coeffs[]) {
    double result = coeffs[0];
    double term = 1.0;
    for (int i = 1; i < n; i++) {
        term *= (x - nodes[i - 1]);
        result += coeffs[i] * term;
    }
    return result;
}

int compare(const void *a, const void *b) {
    double arg1 = *(const double*)a;
    double arg2 = *(const double*)b;
    return (arg1 > arg2) - (arg1 < arg2);
}

int main() {
    int n;
    double a, b;
    
    printf("Введите количество узлов: ");
    scanf("%d", &n);
    printf("Введите границы отрезка [a, b]: ");
    scanf("%lf %lf", &a, &b);
    
    double *eq_nodes = malloc(n * sizeof(double));
    double *ch_nodes = malloc(n * sizeof(double));
    double *eq_values = malloc(n * sizeof(double));
    double *ch_values = malloc(n * sizeof(double));
    
    // Равноотстоящие узлы
    double step = (b - a) / (n - 1);
    for (int i = 0; i < n; i++) {
        eq_nodes[i] = a + i * step;
        eq_values[i] = f(eq_nodes[i]);
    }
    
    // Узлы Чебышева
    for (int i = 0; i < n; i++) {
        double theta = (2 * i + 1) * M_PI / (2 * n);
        ch_nodes[i] = (a + b)/2 + (b - a)/2 * cos(theta);
    }
    qsort(ch_nodes, n, sizeof(double), compare);
    for (int i = 0; i < n; i++) {
        ch_values[i] = f(ch_nodes[i]);
    }
    
    double *eq_coeffs = malloc(n * sizeof(double));
    double *ch_coeffs = malloc(n * sizeof(double));
    computeDividedDifferences(n, eq_nodes, eq_values, eq_coeffs);
    computeDividedDifferences(n, ch_nodes, ch_values, ch_coeffs);
    
    double x_point;
    printf("Введите точку для интерполяции: ");
    scanf("%lf", &x_point);
    
    double eq_result = newtonInterpolate(x_point, n, eq_nodes, eq_coeffs);
    double ch_result = newtonInterpolate(x_point, n, ch_nodes, ch_coeffs);
    double actual_result = f(x_point);
    
    printf("\nРезультат в точке %.3f:\n", x_point);
    printf("  Равноотстоящие узлы: %.6f\n", eq_result);
    printf("  Узлы Чебышева:      %.6f\n", ch_result);
    printf("  Исходная функция:   %.6f\n", actual_result);
    printf("  Отклонения:\n");
    printf("    От равноотстоящих: % .2e\n", eq_result - actual_result);
    printf("    От Чебышева:       % .2e\n\n", ch_result - actual_result);
    
    printf("\nЗначения на отрезке [%.2f, %.2f]:\n", a, b);
    printf("   x      | Равноотстоящие | Чебышева    | Исходная   | Откл.равн | Откл.чеб\n");
    printf("----------|----------------|-------------|------------|-----------|---------\n");
    for (int i = 0; i < 20; i++) {
        double x = a + i * (b - a) / 19;
        double eq_val = newtonInterpolate(x, n, eq_nodes, eq_coeffs);
        double ch_val = newtonInterpolate(x, n, ch_nodes, ch_coeffs);
        double actual_val = f(x);
        
        printf("%-8.3f | %-14.6f | %-11.6f | %-10.6f | %-9.2e | %-9.2e\n", 
              x, eq_val, ch_val, actual_val, 
              eq_val - actual_val, 
              ch_val - actual_val);
    }
    
    free(eq_nodes);
    free(ch_nodes);
    free(eq_values);
    free(ch_values);
    free(eq_coeffs);
    free(ch_coeffs);
    
    return 0;
}
