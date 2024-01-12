#include "utils/vector.cpp"

int main()
{
	int W = 512;
	int H = 512;

	Vector center(0.2, 0.1, 0.);

	std::vector<unsigned char> image(W * H * 3, 0);
#pragma omp parallel for
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < W; j++)
		{

			Vector v(j / (double)W - 0.5, i / (double)H - 0.5, 0.);
			double gaussianVal = exp(-(v - center).norm2() / (2 * sqr(0.2)));

			image[(i * W + j) * 3 + 0] = 127 * gaussianVal; // RED
			image[(i * W + j) * 3 + 1] = 50 * gaussianVal;	// GREEN
			image[(i * W + j) * 3 + 2] = 255 * gaussianVal; // BLUE
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}