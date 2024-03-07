
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
#include "../include/TriangleMesh.hpp"
#include "../include/LightSource.hpp"
#include "../include/Camera.hpp"

#include <iostream>

// Scene createDefaultScene()
// {
// 	Scene scene = Scene();
// 	Sphere lightSource1(Vector(-10, 20, 40), 5, Vector(1, 1, 1), .0, .0, .0, 1e9);
// 	Sphere lightSource2(Vector(10, 20, 40), 5, Vector(1, 1, 1), .0, .0, .0, 1e9);
// 	Sphere lightSource3(Vector(0, 20, 40), 5, Vector(1, 1, 1), .0, .0, .0, 1e9);

// 	scene.addObject(lightSource1);
// 	scene.addObject(lightSource2);
// 	scene.addObject(lightSource3);

// 	Sphere sphere1(Vector(0, 0, 0), 20, Vector(1, 1, 1));					// sphere centrale
// 	Sphere sphere2(Vector(30, 0, 12), 10, Vector(1, 1, 1), 0.0, 0.0, 1.33); // sphere droite
// 	Sphere sphere3(Vector(-30, 0, 12), 10, Vector(1, 1, 1), 1.0);			// sphere gauche

// 	Sphere floor = Sphere(Vector(0, -10000 - 20, 0), 10000, Vector(1, 1, 1));
// 	Sphere ceiling = Sphere(Vector(0, 10000 + 50, 0), 10000, Vector(1, 1, 1));
// 	Sphere wallFront = Sphere(Vector(0, 0, -10000 - 20), 10000, Vector(0, 1, 1));
// 	Sphere wallBack(Vector(0, 0, 10000 + 100), 10000, Vector(1, 0, 1), 0.0);
// 	Sphere wallLeft = Sphere(Vector(-10000 - 50, 0, 0), 10000, Vector(0, 1, 0));
// 	Sphere wallRight = Sphere(Vector(10000 + 50, 0, 0), 10000, Vector(0, 0, 1));

// 	scene.addObject(sphere1);
// 	scene.addObject(sphere2);
// 	scene.addObject(sphere3);
// 	scene.addObject(floor);
// 	scene.addObject(ceiling);
// 	scene.addObject(wallLeft);
// 	scene.addObject(wallRight);
// 	scene.addObject(wallFront);
// 	scene.addObject(wallBack);

// 	return scene;
// }

// Scene myScene()
// {
// 	Scene scene = Scene();
// 	Sphere lightSource(Vector(-10, 20, 40), 10, Vector(1, 1, 1), .0, .0, .0, 3e6);
// 	scene.addObject(lightSource);

// 	Sphere sphere1(Vector(0, 0, 0), 10, Vector(0, 0.6, 0), 0.0, 0.0, 1.5);
// 	Sphere sphere2(Vector(20, 0, -10), 10, Vector(0, 1, 1), 0.0, 0.0, 1.33);
// 	Sphere sphere3(Vector(-30, 0, -10), 10, Vector(1, 0, 1), 1.0);
// 	Sphere floor(Vector(0, -1000, 0), 990, Vector(1, 1, 1), 0.1);
// 	Sphere ceiling(Vector(0, 1000, 0), 940, Vector(1, 0.1, 0.1), 0.0);
// 	Sphere wallBack(Vector(0, 0, 1000), 940, Vector(0.1, 0.1, 1), 0.0);
// 	Sphere wallFront(Vector(0, 0, -1000), 940, Vector(0.1, 0.1, 1), 0.2);
// 	Sphere wallLeft(Vector(-1000, 0, 0), 940, Vector(0.5, 0.1, 1), 0.2);
// 	Sphere wallRight(Vector(1000, 0, 0), 940, Vector(0.1, 0.1, .5), 0.2);

// 	scene.addObject(sphere1);
// 	scene.addObject(sphere2);
// 	scene.addObject(sphere3);
// 	scene.addObject(floor);
// 	scene.addObject(ceiling);
// 	scene.addObject(wallLeft);
// 	scene.addObject(wallRight);
// 	scene.addObject(wallBack);
// 	scene.addObject(wallFront);

// 	return scene;
// }

int main()
{
	bool showProgress = true;
	int s = 512;

	int W = s;
	int H = s;
	Camera camera(Vector(0, 0, 55),  W, H, 80, 1, 55);
	const int nbRays = 20;

	std::vector<unsigned char> image(W * H * 3, 0);

	Scene scene = Scene();
	LightSource lightSource(Vector(-10, 20, 40), 5, 1e9);

	scene.addObject(lightSource);

	Sphere sphere1(Vector(0, 0, 0), 20, Vector(1, 1, 1));					// sphere centrale
	Sphere sphere2(Vector(30, 0, 12), 10, Vector(1, 1, 1), 0.0, 0.0, 1.33); // sphere droite
	Sphere sphere3(Vector(-30, 0, 12), 10, Vector(1, 1, 1), 1.0);			// sphere gauche

	TriangleMesh mesh = TriangleMesh();
	mesh.readOBJ("data/cat.obj");
	Vector barycenter = mesh.getBarycenter();
	mesh.translate(Vector(0, 0, 0) - barycenter);
	mesh.scale(0.05);

	// get bounding box
	std::pair<Vector, Vector> boundingBox = mesh.getBoundingBox();
	std::cout << "Bounding box:" << std::endl;
	std::cout << boundingBox.first.toString() << std::endl;
	std::cout << boundingBox.second.toString() << std::endl;

	Sphere floor = Sphere(Vector(0, -10000 - 20, 0), 10000, Vector(1, 1, 1));
	Sphere ceiling = Sphere(Vector(0, 10000 + 50, 0), 10000, Vector(1, 0, 0));
	Sphere wallFront = Sphere(Vector(0, 0, -10000 - 20), 10000, Vector(0, 1, 1));
	Sphere wallBack(Vector(0, 0, 10000 + 100), 10000, Vector(1, 0, 1), 0.0);
	Sphere wallLeft = Sphere(Vector(-10000 - 50, 0, 0), 10000, Vector(0, 1, 0));
	Sphere wallRight = Sphere(Vector(10000 + 50, 0, 0), 10000, Vector(0, 0, 1));

	// scene.addObject(sphere1);
	// scene.addObject(sphere2);
	// scene.addObject(sphere3);
	scene.addObject(mesh);
	scene.addObject(floor);
	scene.addObject(ceiling);
	scene.addObject(wallLeft);
	scene.addObject(wallRight);
	scene.addObject(wallFront);
	scene.addObject(wallBack);

#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < H; i++)
	{

		for (int j = 0; j < W; j++)
		{
			Vector color(0, 0, 0);
			double di, dj;
			for (int k = 0; k < nbRays; k++)
			{
				Ray ray = camera.launchRay(i, j);
				color += scene.getColor(ray, 1);
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
