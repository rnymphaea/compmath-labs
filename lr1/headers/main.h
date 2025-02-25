#ifndef MAIN_H 
#define MAIN_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double func(double x, double error);
double bisect(double x1, double x2, double length, double error, int * iterations);
double cround(double x, double error);

#endif
