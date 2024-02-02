#pragma once
#include "Vector.hpp"
#include <vector>


double sqr(double x);

Vector computeColor(Vector albedo, Vector L, Vector N, double intensity, double lightVisibility);

std::vector<unsigned char> subSampleImage(std::vector<unsigned char> image, int W, int H, int subSamplingFactor);

double getPercentage(int i, int j, int H, int W);

void boxMuller(double stddev, double &x, double &y);

Vector gammaCorrection(Vector color);
