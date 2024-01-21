#include "../include/Scene.hpp"
#include "../include/Functions.hpp"
#include "../include/Constants.hpp"

#include <cmath>
#include <iostream>

void Scene::add(Sphere s)
{
    spheres.push_back(s);
}

bool Scene::intersect(Ray ray, Vector &P, Vector &N, int &objectIndex)
{
    double t = 1e20;
    bool intersect = false;

    for (size_t i = 0; i < spheres.size(); i++)
    {
        double tTemp = spheres[i].intersectionDistance(ray);
        if (tTemp > 0 && tTemp < t)
        {
            t = tTemp;
            P = ray.origin + t * ray.direction;
            N = (P - spheres[i].center).normalized();
            objectIndex = i;
            intersect = true;
        }
    }
    return intersect;
}

double Scene::isLightVisible(Vector P, Vector L)
{
    Vector lightVector = L - P;
    Ray lightRay(P + lightVector.normalized() * EPSILON, lightVector.normalized());
    double isLightVisible(1);
    for (size_t j = 0; j < spheres.size(); j++)
    {
        double tTemp = spheres[j].intersectionDistance(lightRay);
        if (tTemp > 0 && tTemp < std::sqrt(lightVector.norm2()))
        {
            // si il y a une intersection entre le point et la lumiere avec une autre sphere
            // TODO: tenir compte de la réfraction?
            isLightVisible *= 1 - spheres[j].opacity;
            // break;
        }
    }
    return isLightVisible;
}

Vector Scene::getColor(Vector color, Vector L, Ray ray, double intensity, int depth)
{
    const int maxDepth = 10;
    if (depth > maxDepth)
    {
        return color;
    }
    Vector P, N;
    int intersectIndex = -1;
    bool intersect = this->intersect(ray, P, N, intersectIndex);
    if (intersect)
    {
        Vector lightVector = L - P;

        double isLightVisible = this->isLightVisible(P, L);
        if (isLightVisible < 1e-3)
        {
            return color;
        }
        // si la sphere est reflechissante
        const double reflectance = spheres[intersectIndex].reflectance;
        if (reflectance > 0)
        {
            Vector reflexionVector = ray.direction - 2 * dot(ray.direction, N) * N;
            Ray reflexionRay(P, reflexionVector);
            Vector objectColor = computeColor(isLightVisible * spheres[intersectIndex].albedo, lightVector, N, intensity);
            Vector reflexionColor = getColor(color, L, reflexionRay, intensity, depth + 1);
            return (1 - reflectance) * objectColor + reflectance * reflexionColor;
        }
        else if (spheres[intersectIndex].opacity < 1)
        {
            // entrée ou sortie de la sphère
            bool isEntering = dot(ray.direction, N) < 0;
            double n1, n2;
            if (isEntering)
            {
                n1 = 1.0;
                n2 = spheres[intersectIndex].refractiveIndex;
            }
            else
            {
                n1 = spheres[intersectIndex].refractiveIndex;
                n2 = 1.0;
            }
            double n = n1 / n2;
            double cosThetaI = dot(ray.direction, N);
            Vector tN = std::sqrt(1 - sqr(n) * (1 - sqr(cosThetaI))) * N;
            if (isEntering)
            {
                tN = (-1) * tN;
            }
            Vector tT = n * ray.direction;

            Ray refractionRay = Ray(P + (tT + tN).normalized() * EPSILON, (tN + tT).normalized());
            Vector refractionColor = getColor(color, L, refractionRay, intensity, depth + 1);
            return refractionColor;
        }
        else
        {
            Vector objectColor = computeColor(isLightVisible * spheres[intersectIndex].albedo, lightVector, N, intensity);
            return color + objectColor;
        }
    }
    return color;
}
