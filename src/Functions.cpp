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

char *imageName(const char *name, const char *ext, int resolution, int nbRays, int renderTime)
{
	char resolutionString[20];
	char raysString[20];
	char timeString[20];

	sprintf(resolutionString, "%dp", resolution);
	sprintf(raysString, "%dr", nbRays);

	if (renderTime < 10000)
	{
		sprintf(timeString, "%dms", renderTime);
	}
	else
	{
		int seconds = renderTime / 1000;
		sprintf(timeString, "%ds", seconds);
	}

	char *result = new char[strlen(name) + strlen(resolutionString) + strlen(raysString) + strlen(timeString) + 4];
	strcpy(result, name);
	strcat(result, "_");
	strcat(result, resolutionString);
	strcat(result, "_");
	strcat(result, raysString);
	strcat(result, "_");
	strcat(result, timeString);
	strcat(result, ".");
	strcat(result, ext);

	return result;
}
