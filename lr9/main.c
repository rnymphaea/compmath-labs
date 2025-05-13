#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_ITER 30

double f(double x) {
    return cos(x*x);
}

double rectangle(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        double x = a + (i + 0.5) * h;
        sum += f(x);
    }
    return h * sum;
}

double trapezoid(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.5 * (f(a) + f(b));
    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        sum += f(x);
    }
    return h * sum;
}

double simpson(double a, double b, int n) {
    if (n % 2 != 0) n++;
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        sum += (i % 2 == 0) ? 2.0 * f(x) : 4.0 * f(x);
    }
    return h * sum / 3.0;
}

void compute_integral(double a, double b, double eps, double (*method)(double, double, int), int method_type, const char* name) {
    int n = (method_type == 2) ? 2 : 1;
    double I_prev, I_curr, delta, refined_I, error;
    int best_n = n;
    
    printf("\nMethod: %s\n", name);
    printf("--------------------------------\n");
    
    do {
        I_prev = method(a, b, n);
        I_curr = method(a, b, 2 * n);
        delta = fabs(I_curr - I_prev);
        
        switch (method_type) {
            case 0: // Прямоугольники
                refined_I = I_curr + delta / 3.0;
                error = delta / 3.0;
                break;
            case 1: // Трапеции
                refined_I = I_curr - delta / 3.0;
                error = delta / 3.0;
                break;
            case 2: // Симпсон
                refined_I = I_curr - delta / 15.0;
                error = delta / 15.0;
                break;
        }
        
        printf("n = %4d | Integral = %.10f | Error = %.3e\n", n, refined_I, error);
        
        best_n = n;
        n *= 2;
    } while (error > eps && n < (1 << MAX_ITER));
    
    printf("\nResult: %.10f\nRequired n: %d\n", refined_I, best_n);
}

int main() {
    double a = 0.0;
    double b = 1.0;
    double eps;

    printf("Введите требуемую точность (например, 1e-6): ");
    if (scanf("%lf", &eps) != 1 || eps <= 0) {
        printf("Ошибка ввода! Допустимы только положительные числа\n");
        return 1;
    }

    compute_integral(a, b, eps, rectangle, 0, "Rectangle");
    compute_integral(a, b, eps, trapezoid, 1, "Trapezoid");
    compute_integral(a, b, eps, simpson, 2, "Simpson");

    return 0;
}
