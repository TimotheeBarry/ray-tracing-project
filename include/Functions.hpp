#pragma once
#include "Vector.hpp"
#include <vector>
#include <random>


double sqr(double x);

void boxMuller(double stddev, double &x, double &y);

Vector gammaCorrection(Vector color);
