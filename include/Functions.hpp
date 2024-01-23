#pragma once
#include "Vector.hpp"
#include <vector>


double sqr(double x);

Vector computeColor(Vector albedo, Vector L, Vector N, double intensity, double lightVisibility);

std::vector<unsigned char> subSampleImage(std::vector<unsigned char> image, int W, int H, int subSamplingFactor);

double getPercentage(int i, int j, int H, int W);

Vector gammaCorrection(Vector color);
