
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
#include "../include/BoundingBox.hpp"

#include <iostream>
#include <chrono>

int main()
{

	bool showProgress = true;
	int s = 512;

	int W = s;
	int H = s;
	Camera camera(Vector(0, 0, 55), W, H, 60, 1, 55);
	const int nbRays = 25;

	std::vector<unsigned char> image(W * H * 3, 0);

	Scene scene = Scene();
	LightSource lightSource(Vector(-10, 20, 40), 5, 4e9);

	scene.addObject(lightSource);

	Diffuse diffuse;		 // materiau diffus (lambertien)
	Transparent glass(1.52); // verre
	Mirror mirror(1);		 // miroir
	BlinnPhong blinnPhong(1000, 0.1);

	Sphere floor = Sphere(Vector(0, -10000 - 20, 0), 10000, Vector(1, 1, 1), &diffuse);
	Sphere ceiling = Sphere(Vector(0, 10000 + 50, 0), 10000, Vector(1, 0, 0), &diffuse);
	Sphere wallFront = Sphere(Vector(0, 0, -10000 - 50), 10000, Vector(0, 1, 1), &diffuse);
	Sphere wallBack(Vector(0, 0, 10000 + 100), 10000, Vector(1, 0, 1), &diffuse);
	Sphere wallLeft = Sphere(Vector(-10000 - 50, 0, 0), 10000, Vector(0, 1, 0), &diffuse);
	Sphere wallRight = Sphere(Vector(10000 + 50, 0, 0), 10000, Vector(0, 0, 1), &diffuse);

	Sphere sphere1(Vector(0, 20, 0), 5, Vector(1, .0, 0), &blinnPhong); // sphere centrale
	Sphere sphere2(Vector(20, 15, -10), 5, Vector(1, .5, .3), &mirror); // sphere droite
	Sphere sphere3(Vector(-10, 12, 25), 5, Vector(1, 1, 1), &glass);	// sphere gauche

	TriangleMesh mesh = TriangleMesh("data/cat.obj", "data/cat_diff.png", &blinnPhong);

	Vector barycenter = mesh.getBarycenter();
	mesh.translate(Vector(0, 0, 0) - barycenter);
	mesh.scale(0.8);
	mesh.translate(Vector(0, floor.center[1] + floor.radius - mesh.bbox.min[1], 0));
	mesh.rotate(-PI / 4, Vector(0, 1, 0));
	mesh.initBVH();

	scene.addObject(sphere1);
	scene.addObject(sphere2);
	scene.addObject(sphere3);
	scene.addObject(floor);
	scene.addObject(ceiling);
	scene.addObject(wallLeft);
	scene.addObject(wallRight);
	scene.addObject(wallFront);
	scene.addObject(wallBack);
	scene.addObject(mesh);

	auto start = std::chrono::high_resolution_clock::now();
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
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

	return 0;
}
