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
	double eps;
	double ye;
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
		printf("error: invalid delta!\n");
		exit(1);
	}
	if (x > 0) {
		return error * (int)(x / error + 0.5);
	}

	return error * (int)(x / error - 0.5);
}

int main() {
	double left = 0;
	double right = 0.5;
	double length;
	double error = 0.000001;
	int iterations;
	for (int step_error = 0; step_error < 6; step_error++) {
		length = 0.1;
		printf("\t\t\t\terror = %.6f\n", error);
		printf(ROW);
		printf("|\tlength\t\t|\tsolution (x)\t|\ty(x)\t\t|   iterations\t|\n");
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
	return 0;
}

