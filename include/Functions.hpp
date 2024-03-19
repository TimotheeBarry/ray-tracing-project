#pragma once
#include "Vector.hpp"
#include <vector>
#include <random>
#include <cstring>


double sqr(double x);

void boxMuller(double stddev, double &x, double &y);

Vector gammaCorrection(Vector color);

char *imageName(const char *name, const char *ext,int resolution, int nbRays, int renderTime);