
#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../stb/stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "../stb/stb_image.h"

#include "../include/Vector.hpp"
#include "../include/Ray.hpp"
#include "../include/Sphere.hpp"
#include "../include/Scene.hpp"
#include "../include/Functions.hpp"
#include "../include/Constants.hpp"

#include <iostream>

Scene createDefaultScene()
{
	Scene scene = Scene();
	Sphere lightSource1(Vector(-10, 20, 40), 5, Vector(1, .0, .0), .0, .0, .0, 1e9);
	Sphere lightSource2(Vector(10, 20, 40), 5, Vector(.0, 1, .0), .0, .0, .0, 1e9);
	Sphere lightSource3(Vector(0, 20, 40), 5, Vector(.0, .0, 1), .0, .0, .0, 1e9);
	
	scene.addSphere(lightSource1);
	scene.addSphere(lightSource2);
	scene.addSphere(lightSource3);

	Sphere sphere1(Vector(0, 0, 0), 20, Vector(1, 1, 1)); // sphere centrale
	Sphere sphere2(Vector(30, 0, 12), 10, Vector(1, 1, 1), 0.0, 0.0, 1.33); // sphere droite
	Sphere sphere3(Vector(-30, 0, 12), 10, Vector(1, 1, 1), 1.0); // sphere gauche

	Sphere floor = Sphere(Vector(0, -10000 - 20, 0), 10000, Vector(1, 1, 1));
	Sphere ceiling = Sphere(Vector(0, 10000 + 50, 0), 10000, Vector(1, 1, 1));
	Sphere wallFront = Sphere(Vector(0, 0, -10000 - 20), 10000, Vector(0, 1, 1));
	Sphere wallLeft = Sphere(Vector(-10000 - 50, 0, 0), 10000, Vector(0, 1, 0));
	Sphere wallRight = Sphere(Vector(10000 + 50, 0, 0), 10000, Vector(0, 0, 1));

	scene.addSphere(sphere1);
	scene.addSphere(sphere2);
	scene.addSphere(sphere3);
	scene.addSphere(floor);
	scene.addSphere(ceiling);
	scene.addSphere(wallLeft);
	scene.addSphere(wallRight);
	scene.addSphere(wallFront);

	return scene;
}

Scene myScene()
{
	Scene scene = Scene();
	Sphere lightSource(Vector(-10, 20, 40), 10, Vector(1, 1, 1), .0, .0, .0, 3e6);
	scene.addSphere(lightSource);

	Sphere sphere1(Vector(0, 0, 0), 10, Vector(0, 0.6, 0), 0.0, 0.0, 1.5);
	Sphere sphere2(Vector(20, 0, -10), 10, Vector(0, 1, 1), 0.0, 0.0, 1.33);
	Sphere sphere3(Vector(-30, 0, -10), 10, Vector(1, 0, 1), 1.0);
	Sphere floor(Vector(0, -1000, 0), 990, Vector(1, 1, 1), 0.1);
	Sphere ceiling(Vector(0, 1000, 0), 940, Vector(1, 0.1, 0.1), 0.0);
	Sphere wallBack(Vector(0, 0, 1000), 940, Vector(0.1, 0.1, 1), 0.0);
	Sphere wallFront(Vector(0, 0, -1000), 940, Vector(0.1, 0.1, 1), 0.2);
	Sphere wallLeft(Vector(-1000, 0, 0), 940, Vector(0.5, 0.1, 1), 0.2);
	Sphere wallRight(Vector(1000, 0, 0), 940, Vector(0.1, 0.1, .5), 0.2);

	scene.addSphere(sphere1);
	scene.addSphere(sphere2);
	scene.addSphere(sphere3);
	scene.addSphere(floor);
	scene.addSphere(ceiling);
	scene.addSphere(wallLeft);
	scene.addSphere(wallRight);
	scene.addSphere(wallBack);
	scene.addSphere(wallFront);

	return scene;
}

int main()
{
	bool showProgress = true;
	int s = 512;
	int W = s;
	int H = s;

	double alpha = 80 * (PI) / 180;
	const int nbRays = 25;

	std::vector<unsigned char> image(W * H * 3, 0);

	Vector O(0, 0, 55); // camera origin

	Scene scene = createDefaultScene();
	// Scene scene = myScene();

#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < H; i++)
	{

		for (int j = 0; j < W; j++)
		{
			Vector color(0, 0, 0);
			double di, dj;
			for (int k = 0; k < nbRays; k++)
			{

				double di, dj;
				boxMuller(0.2, di, dj);
				Vector vector(j - W / 2 + 0.5 + dj, -i + H / 2 - 0.5 + di, -W / (2 * std::tan(alpha / 2)));
				Ray ray(O, vector.normalized());
				color += scene.getColor(ray, 5);
			}
			color /= nbRays;

			color.clip(0, 255);
			color = gammaCorrection(color);

			image[(i * W + j) * 3 + 0] = color[0]; // RED
			image[(i * W + j) * 3 + 1] = color[1]; // GREEN
			image[(i * W + j) * 3 + 2] = color[2]; // BLUE
		}
	}

	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}
