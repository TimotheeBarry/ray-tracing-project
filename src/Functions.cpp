#include "../include/Functions.hpp"
#include "../include/Constants.hpp"
#include <cmath>

double sqr(double x)
{
    return x * x;
}

Vector gammaCorrection(Vector color)
{
    return pow(color / 255, 1 / GAMMA)*255;
}

Vector computeColor(Vector albedo, Vector L, Vector N, double intensity, double lightVisibility)
{
    return lightVisibility * albedo * intensity * (std::max(0., dot(L.normalized(), N)) / (4 * sqr(PI) * L.norm2()));
}

std::vector<unsigned char> subSampleImage(std::vector<unsigned char> image, int W, int H, int subSamplingFactor)
{
	if (subSamplingFactor == 1)
	{
		return image;
	}
	int W_bis = W / subSamplingFactor;
	int H_bis = H / subSamplingFactor;

	std::vector<unsigned char> subSampledImage(W_bis * H_bis * 3, 0);

	for (int i = 0; i < H; i += subSamplingFactor)
	{
		for (int j = 0; j < W; j += subSamplingFactor)
		{
			int index = (i * W + j) * 3;
			int subSampledIndex = (i / subSamplingFactor * W_bis + j / subSamplingFactor) * 3;
			for (int k = 0; k < 3; k++)
			{
				int sum = 0;
				for (int l = 0; l < subSamplingFactor; l++)
				{
					for (int m = 0; m < subSamplingFactor; m++)
					{
						sum += image[index + (l * W + m) * 3 + k];
					}
				}
				subSampledImage[subSampledIndex + k] = sum / (subSamplingFactor * subSamplingFactor);
			}
		}
	}

	return subSampledImage;
}

//fonction qui donne le pourcentage de remplissage d'un carré de taille H*W pour i et j
double getPercentage(int i, int j, int H, int W)
{
    return (double)(i * W + j) / (H * W);
}
