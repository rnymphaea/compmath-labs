// var 14: f(x) = cos(x^2) - 10x

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double func(double x, double error);
double chord(double x1, double x2, double eps, double error, int * iterations);
double cround(double x, double error);


#define ROW "+-----------------------+-----------------------+-----------------------+---------------+\n"

double func(double x, double error) {
	double result = cos(x*x) - 10*x;
	if (error != 0.000001) {
		return cround(result, error);
	}
	return result;
}

double chord(double x1, double x2, double eps, double error, int * iterations) {
	double f1 = func(x1, error);
	if (f1 == 0) {
		return x1;
	}

	double f2 = func(x2, error);
	if (f2 == 0) {
		return x2;
	}

	double c, fc;
	do {
		c = x1 - ( (x2 - x1) * f1 / (f2 - f1) );
		fc = func(c, error);
		if (fc == 0) {
			return c;
		}

		if (fc * f1 < 0) {
			x2 = c;
			f2 = fc;
		}
		else {
			x1 = c;
			f1 = fc;
		}
		(*iterations)++;
	} while (fabs(fc) >= eps);
	
	return c;
}

double cround(double x, double error) {
	if (error < 1E-9) {
		printf("invalid error!\n");
		exit(1);
	}
	if (x > 0) {
		return error * (int)(x / error + 0.5);
	}

	return error * (int)(x / error - 0.5);
}

int main() {
	double left = 0; // левая граница начального отрезка 
	double right = 1; // правая граница начального отрезка 
	double eps; // радиус отрезка, в который должен попасть func(root)

	double error = 0.000001; // значение ошибки для округления
	int iterations; // количество итераций

	for (int step_error = 0; step_error < 6; step_error++) {
		eps = 0.1;

		printf("\t\t\t\terror = %.6f\n", error);
		printf(ROW);
		printf("|\t eps\t\t|\tsolution (root)\t|\t y(root)\t|   iterations\t|\n");
		printf(ROW);

		for (int step = 0; step < 6; step++) {
			iterations = 0;
			double result = chord(left, right, eps, error, &iterations);

			printf("|\t%.6f\t|\t%.6f\t|\t", eps, result);
			printf("%.6f\t|\t",  func(result, error));
			printf("%d\t|\n", iterations);
			
			eps *= 0.1;
		}

		error /= 0.1;

		printf(ROW);
		printf("\n\n");
	}
	
	printf("error - значение ошибки для округления\n");
        printf("eps - радиус отрезка, в который должен попасть func(root)\n");
        printf("solution (root) - корень\n");
        printf("y(root) - значение функции в точке root\n");
        printf("iterations - количество итераций\n");

	return 0;
}

