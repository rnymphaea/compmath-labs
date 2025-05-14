#include <stdio.h>
#include <math.h>


double f(double x) {
    return sin(x) / (1.0 + x);
}

double gauss_integral() {
    const double X[] = {
        -0.96028986, -0.79666648, -0.52533242, -0.18343464,
         0.18343464,  0.52533242,  0.79666648,  0.96028986
    };
    
    const double A[] = {
        0.10122854, 0.22238103, 0.31370664, 0.36268378,
        0.36268378, 0.31370664, 0.22238103, 0.10122854
    };

    double sum = 0.0;
    const double a = 0.0, b = 1.0; // Пределы интегрирования
    
    // Преобразование узлов в интервал [a, b] и вычисление суммы
    for (int i = 0; i < 8; i++) {
        double t = 0.5 * (X[i] * (b - a) + (a + b));
        sum += A[i] * f(t);
    }
    
    return 0.5 * (b - a) * sum; // Масштабирование результата
}

int main() {
    double result = gauss_integral();
    
    printf("Result:\n");
    printf("∫sin(x)/(1+x)dx ≈ %.15f\n", result);
    
    return 0;
}
