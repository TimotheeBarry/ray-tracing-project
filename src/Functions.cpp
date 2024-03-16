#include "../include/Functions.hpp"
#include "../include/Constants.hpp"
#include <cmath>

double sqr(double x)
{
	return x * x;
}

Vector gammaCorrection(Vector color)
{
	return pow(color / 255, 1 / GAMMA) * 255;
}


void boxMuller(double stddev, double &x, double &y)
{
	double r1 = uniform(gen);
	double r2 = uniform(gen);

	x = sqrt(-2 * log(r1)) * cos(2 * PI * r2) * stddev;
	y = sqrt(-2 * log(r1)) * sin(2 * PI * r2) * stddev;
}