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

double Scene::lightVisibility(Vector P, Vector lightSource)
{
    Vector L = lightSource - P;
    Ray lightRay(P + L.normalized() * EPSILON, L.normalized());
    double lightVisibility(1);
    for (size_t j = 0; j < spheres.size(); j++)
    {
        double tTemp = spheres[j].intersectionDistance(lightRay);
        if (tTemp > 0 && tTemp < std::sqrt(L.norm2()))
        {
            // si il y a une intersection entre le point et la lumiere avec une autre sphere
            // TODO: tenir compte de la réfraction?
            lightVisibility *= 1 - spheres[j].opacity;
            // break;
        }
    }
    return lightVisibility;
}

Vector Scene::getColor(Vector color, Vector lightSource, Ray ray, double intensity, int depth)
{
    const int maxDepth = 8;
    if (depth > maxDepth)
    {
        return color;
    }
    Vector P, N;
    int intersectIndex = -1;
    bool intersect = this->intersect(ray, P, N, intersectIndex);
    if (intersect)
    {
        Sphere sphere = spheres[intersectIndex];
        Vector L = lightSource - P; // light vector from point P to light source

        double lightVisibility = this->lightVisibility(P, lightSource);
        // si la sphere est reflechissante
        const double reflectance = sphere.reflectance;
        if (reflectance > 0)
        {
            // couleur diffusée
            Vector diffusedColor = computeColor(sphere.albedo, L, N, intensity, lightVisibility);
            // couleur réfléchie
            Vector reflexionVector = ray.direction - 2 * dot(ray.direction, N) * N;
            Ray reflectedRay(P + N * EPSILON, reflexionVector.normalized());
            Vector reflectedColor = getColor(color, lightSource, reflectedRay, intensity, depth + 1);

            return (1 - reflectance) * diffusedColor + reflectance * reflectedColor;
        }
        else if (sphere.opacity < 1)
        {
            // entrée ou sortie de la sphère
            bool isEntering = dot(ray.direction, N) < 0;
            double n1, n2;
            if (isEntering)
            {
                n1 = 1.0;
                n2 = sphere.refractiveIndex;
                N = (-1) * N;
            }
            else
            {
                n1 = sphere.refractiveIndex;
                n2 = 1.0;
            }
            double n = n1 / n2;
            double cosThetaI = dot(ray.direction, N);
            Vector tN = std::sqrt(1 - sqr(n) * (1 - sqr(cosThetaI))) * N;

            Vector tT = n * ray.direction;

            Ray refractionRay = Ray(P + N * EPSILON, (tN + tT).normalized());
            Vector refractionColor = getColor(color, lightSource, refractionRay, intensity, depth + 1);
            double lightTransmissionCoef = isEntering ? sphere.opacity : 1;
            return lightTransmissionCoef * refractionColor;
        }
        else
        {
            Vector diffusedColor = computeColor(sphere.albedo, L, N, intensity, lightVisibility);
            return diffusedColor;
        }
    }
    return color;
}
