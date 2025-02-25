// var 14: f(x) = cos(x^2) - 10x

#include "./headers/main.h"

#define ROW "+-----------------------+-----------------------+-----------------------+---------------+\n"

double func(double x, double error) {
	double result = cos(x*x) - 10*x;
	if (error != 0.000001) {
		return cround(result, error);
	}
	return result;
}

double bisect(double x1, double x2, double length, double error, int * iterations) {
	double eps; // корень 
	double ye; // значение функции в точке eps
	while (x2 - x1 >= 2*length) {
		eps = (x1 + x2) / 2;
		ye = func(eps, error);
		if (ye == 0) {
			return eps;
		}
		if (func(x1, error) * ye < 0) {
			x2 = eps;	
		} else {
			x1 = eps;
		}
		(*iterations)++;
	}
	return eps;
}

double cround(double x, double error) {
	if (error <= 0.000001) {
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
	double right = 0.5; // правая граница начального отрезка 
	double length; // длина отрезка, в который должен попасть eps, чтобы быть корнем

	double error = 0.000001; // значение ошибки для округления
	int iterations; // количество итераций

	for (int step_error = 0; step_error < 6; step_error++) {
		length = 0.1;

		printf("\t\t\t\terror = %.6f\n", error);
		printf(ROW);
		printf("|\t length\t\t|\tsolution (eps)\t|\t y(eps)\t\t|   iterations\t|\n");
		printf(ROW);

		for (int step = 0; step < 6; step++) {
			iterations = 0;
			double result = bisect(left, right, length, error, &iterations);

			printf("|\t%.6f\t|\t%.6f\t|\t", length, result);
			printf("%.6f\t|\t",  func(result, error));
			printf("%d\t|\n", iterations);
			
			length *= 0.1;
		}

		error /= 0.1;

		printf(ROW);
		printf("\n\n");
	}
	
	printf("error - значение ошибки для округления\n");
        printf("length - длина отрезка, в который должен попасть eps, чтобы быть корнем\n");
        printf("solution (eps) - корень\n");
        printf("y(eps) - значение функции в точке eps\n");
        printf("iterations - количество итераций\n");

	return 0;
}

