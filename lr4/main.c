// var 14: f(x) = cos(x^2) - 10x

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

double func(double x, double error);
double fi(double x, double lambda, double error);
double derivative_fi(double x, double lambda, double error);
double iter(double x, double eps, double lambda, double error, int * iterations);
double cround(double x, double error);


#define ROW "+-----------------------+-----------------------+-----------------------+---------------+\n"

double func(double x, double error) {
	double result = cos(x*x) - 10*x;
	if (error != 0.000001) {
		return cround(result, error);
	}
	return result;
}

double fi(double x, double lambda, double error) {
	double result = x - lambda * (cos(x*x) - 10*x);
	if (error != 0.000001) {
		return cround(result, error);
	}
	return result;
}

double derivative_fi(double x, double lambda, double error) {
	double result = 2 * lambda * x * sin(x*x) + 10*lambda + 1;
	if (error != 0.000001) {
		return cround(result, error);
	}
	return result;
}

double iter(double x, double eps, double lambda, double error, int * iterations) {
	double x1 = fi(x, lambda, error);
	double x2 = fi(x1, lambda, error);
	(*iterations) = 2;
	
	while ( (x1 - x2) * (x1 - x2) >= fabs((2 * x1 - x - x2) * eps) ) {
		x = x1;
		x1 = x2;
		x2 = fi(x1, lambda, error);
		(*iterations)++;
	}
	
	return x2;
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
	double start = 0.1;
	double eps = 0.000001;
	double error = 0.000001;
	double lambda = -0.098778;

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
				double result = iter(start, eps, lambda, error, &iterations);

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
	}

	else {
		int opt = getopt(argc, argv, "s:e:E:l:h");
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
				case 'l':
					lambda = atof(optarg);
					break;
				case 'h':
					printf("'-s' - стартовое значение\n'-e' - значение ошибки для округления\n'-E' - значение eps\n'-l' - значение коэффициента при функции для fi(x)");
					return 0;
						
			}
			opt = getopt(argc, argv, "s:e:E:l:h");
		}
		
		printf("\t\t\t\terror = %.6f\n", error);
		printf(ROW);
		printf("|\t eps\t\t|\tsolution (root)\t|\t y(root)\t|   iterations\t|\n");
		printf(ROW);

		double result = iter(start, eps, lambda, error, &iterations);
		printf("|\t%.6f\t|\t%.6f\t|\t", eps, result);
                printf("%.6f\t|\t",  func(result, error));
                printf("%d\t|\n", iterations);
		printf(ROW);
		
	}
	return 0;
}

