
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
	// Materiaux principaux
	Diffuse diffuse;
	Transparent glass(1.52);
	Transparent water(1.33);
	Transparent diamond(2.42);
	Mirror mirror(1);
	Mirror semiMirror(.1);
	Mirror semiMirror2(.5);
	BlinnPhong blinnPhong_001(1000, 0.01);
	BlinnPhong blinnPhong_0025(1000, 0.025);
	BlinnPhong blinnPhong_005(1000, 0.05);
	BlinnPhong blinnPhong_010(1000, 0.1);

	// Objets principaux
	LightSource lightSource(Vector(-10, 20, 40), 5, 4e9); // default -10 20 40
	Sphere floor = Sphere(Vector(0, -10000 - 20, 0), 10000, Vector(1, 1, 1), &diffuse);
	Sphere ceiling = Sphere(Vector(0, 10000 + 50, 0), 10000, Vector(1, 0.01, 0.01), &diffuse);
	Sphere wallFront = Sphere(Vector(0, 0, -10000 - 50), 10000, Vector(0.01, 1, 1), &diffuse);
	Sphere wallBack(Vector(0, 0, 10000 + 100), 10000, Vector(1, 0.01, 1), &diffuse);
	Sphere wallLeft = Sphere(Vector(-10000 - 50, 0, 0), 10000, Vector(0.01, 1, 0.01), &diffuse);
	Sphere wallRight = Sphere(Vector(10000 + 50, 0, 0), 10000, Vector(0.01, 0.01, 1), &diffuse);

	// Camera
	int res = 256;
	int nbRays = 200;
	int fov = 60;
	int aperture = 0;
	int focalLength = 55;
	Camera cam(Vector(0, 0, 55), res, res, nbRays, fov, aperture, focalLength);

	// Scene - Murs et sol
	Scene scene = Scene();
	scene.addObject(floor);
	scene.addObject(ceiling);
	scene.addObject(wallLeft);
	scene.addObject(wallRight);
	scene.addObject(wallFront);
	scene.addObject(wallBack);
	scene.addObject(lightSource);

	// Scene - Objets principaux

	Sphere sphere1 = Sphere(Vector(-20, -10, -20), 10, Vector(1, 0.05, 0.05), &diffuse);
	Sphere sphere2 = Sphere(Vector(10, 10, 20), 10, Vector(0.1, 1, .3), &diffuse);
	Sphere sphere3 = Sphere(Vector(0, 0, 0), 20, Vector(0.1, .3, 1), &diffuse);
	// TriangleMesh mesh1 = TriangleMesh("data/cat.obj", "data/cat_diff.png", &diffuse);
	// TriangleMesh mesh2 = TriangleMesh("data/cat.obj", "data/cat_diff.png", &diffuse);
	// TriangleMesh mesh3 = TriangleMesh("data/cat.obj", "data/cat_diff.png", &glass);

	// Vector barycenter1 = mesh1.getBarycenter();
	// mesh1.translate(Vector(-16, 0, 0) - barycenter1);
	// mesh1.scale(.4);
	// mesh1.rotate(PI / 4, Vector(0, 0, 1));
	// mesh1.translate(Vector(0, floor.center[1] + floor.radius - mesh1.bbox.min[1], 0));
	// mesh1.rotate(-PI / 4, Vector(0, 1, 0));
	// mesh1.initBVH();

	// Vector barycenter2 = mesh2.getBarycenter();
	// mesh2.translate(Vector(16, 0, 0) - barycenter2);
	// mesh2.scale(.8);
	// mesh2.translate(Vector(0, floor.center[1] + floor.radius - mesh2.bbox.min[1], 0));
	// mesh2.rotate(-3 * PI / 4, Vector(0, 1, 0));
	// mesh2.initBVH();

	// Vector barycenter3 = mesh3.getBarycenter();
	// mesh3.translate(Vector(16, 0, 0) - barycenter3);
	// mesh3.scale(.7);
	// mesh3.rotate(-PI / 3, Vector(0, 1, 0));
	// mesh3.translate(Vector(0, floor.center[1] + floor.radius - mesh3.bbox.min[1], 0));

	// mesh3.initBVH();

	// Scene - Ajout des objets à la scène
	// scene.addObject(mesh1);
	// scene.addObject(mesh2);
	// scene.addObject(mesh3);
	// scene.addObject(sphere1);
	// scene.addObject(sphere2);
	scene.addObject(sphere3);

	// Génération de l'image
	std::vector<unsigned char> image(cam.width * cam.height * 3, 0);
	auto start = std::chrono::high_resolution_clock::now();
#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < cam.height; i++)
	{
		for (int j = 0; j < cam.width; j++)
		{
			Vector color(0, 0, 0);
			double di, dj;
			for (int k = 0; k < nbRays; k++)
			{
				Ray ray = cam.launchRay(i, j);
				color += scene.getColor(ray, 10);
			}
			color /= nbRays;

			color.clip(0, 255);
			color = gammaCorrection(color);

			image[(i * cam.width + j) * 3 + 0] = color[0]; // RED
			image[(i * cam.width + j) * 3 + 1] = color[1]; // GREEN
			image[(i * cam.width + j) * 3 + 2] = color[2]; // BLUE
		}
	}

	auto end = std::chrono::high_resolution_clock::now();
	int duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000. << "s" << std::endl;
	auto resultName = imageName("generated/image", "png", res, nbRays, duration);
	stbi_write_png(resultName, cam.width, cam.height, 3, &image[0], 0);

	return 0;
}
