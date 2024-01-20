#pragma once
#include "Vector.hpp"

double sqr(double x);

Vector computeColor(Vector albedo, Vector lightVector, Vector N, double intensity);

Vector gammaCorrection(Vector color);
