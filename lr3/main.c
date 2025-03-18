// var 14: f(x) = cos(x^2) - 10x

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

double func(double x, double error);
double derivative(double x, double error);
double newton(double x, double eps, double error, int * iterations);
double cround(double x, double error);


#define ROW "+-----------------------+-----------------------+-----------------------+---------------+\n"

double func(double x, double error) {
	double result = cos(x*x) - 10*x;
	if (error != 0.000001) {
		return cround(result, error);
	}
	return result;
}

double derivative(double x, double error) {
	double result = -2 * x * sin(x*x) - 10;
	if (error != 0.000001) {
		return cround(result, error);
	}
	return result;
}

double newton(double x, double eps, double error, int * iterations) {
	double y, d;
	
	do {
		y = func(x, error);
		if (y == 0) {
			return x;	
		}

		d = derivative(x, error);
		if (d == 0) {
			printf("derivative = 0\n");
			exit(1);
		}
		
		x -= y / d;

		(*iterations)++;
	} while (fabs(y / d) > eps);

	return x;
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

int main(int argc, char * argv[]) {
	double start = 0.25;
	double eps = 0.000001;
	double error = 0.000001;

	int iterations = 0; // количество итераций
	
	if (argc == 1) {
		for (int step_error = 0; step_error < 6; step_error++) {
			eps = 0.1;

			printf("\t\t\t\terror = %.6f\n", error);
			printf(ROW);
			printf("|\t eps\t\t|\tsolution (root)\t|\t y(root)\t|   iterations\t|\n");
			printf(ROW);

			for (int step = 0; step < 6; step++) {
				iterations = 0;
				double result = newton(start, eps, error, &iterations);

				printf("|\t%.6f\t|\t%.8f\t|\t", eps, result);
				printf("%.8f\t|\t",  func(result, error));
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
	}

	else {
		int opt = getopt(argc, argv, "s:e:E:h");
		while (opt != -1) {
			switch (opt) {
				case 's':
					start = atof(optarg);
					break;
				case 'e':
					error = atof(optarg);
					break;
				case 'E':
					eps = atof(optarg);
					break;
				case 'h':
					printf("'-s' - стартовое значение\n'-e' - значение ошибки для округления\n'-E' - значение eps\n");
					return 0;
						
			}
			opt = getopt(argc, argv, "s:e:E:h");
		}
		
		printf("\t\t\t\terror = %.6f\n", error);
		printf(ROW);
		printf("|\t eps\t\t|\tsolution (root)\t|\t y(root)\t|   iterations\t|\n");
		printf(ROW);

		double result = newton(start, eps, error, &iterations);
		printf("|\t%.6f\t|\t%.6f\t|\t", eps, result);
                printf("%.6f\t|\t",  func(result, error));
                printf("%d\t|\n", iterations);
		printf(ROW);
		
	}
	return 0;
}

