#include "../include/Scene.hpp"
#include "../include/Functions.hpp"
#include "../include/Constants.hpp"
#include <cmath>
#include <iostream>

void Scene::addSphere(Sphere s)
{
    spheres.push_back(s);
}

void Scene::addLight(Sphere s)
{
    lights.push_back(s);
}

bool Scene::intersect(Ray &ray, Vector &P, Vector &N, int &objectIndex, double &t)
{
    t = 1e20;
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

Vector Scene::getColor(Ray &ray, int depth)
{
    if (depth < 0)
    {
        return Vector(0, 0, 0);
    }

    Vector P, N;
    int intersectIndex = -1;
    double t;
    bool intersect = this->intersect(ray, P, N, intersectIndex, t);
    if (intersect)
    {

        Sphere &sphere = spheres[intersectIndex];
        const double reflectance = sphere.reflectance;
        const double opacity = sphere.opacity;
        double R(-1), T(-1); // coeffecient de transmission

        // Si la sphere est transparente
        if (opacity < 1)
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
            double k0 = sqr((n1 - n2) / (n1 + n2));

            R = k0 + (1 - k0) * std::pow(1 - std::abs(dot(ray.direction, N)), 5);
            T = 1 - R;

            double n = n1 / n2;
            double cosThetaI = dot(ray.direction, N);
            Vector tN = std::sqrt(1 - sqr(n) * (1 - sqr(cosThetaI))) * N;

            Vector tT = n * ray.direction;
            // on choisit aléatoirement entre réflexion et transmission avec proba qui dépend des coefficients R et T
            if (uniform(gen) < T)
            {
                // transmission
                Ray ray = Ray(P + N * EPSILON, (tN + tT).normalized());
                return getColor(ray, depth - 1);
            }
            else
            {
                // réflexion
                Vector reflexionVector = ray.direction - 2 * dot(ray.direction, N) * N;
                Ray reflectedRay(P + N * EPSILON, reflexionVector.normalized());
                return getColor(reflectedRay, depth - 1);
            }
        }
        else
        {
            Vector indirectColor(0, 0, 0), diffusedColor(0, 0, 0);
            // on calcule la contribution directe et indirecte que si la sphère n'est pas un miroir pur pour éviter les calculs inutiles
            if (reflectance < 1)
            {
                // contribution indirecte
                Vector randomVector = N.generateRandomCosineVector();
                Ray randomRay = Ray(P + N * EPSILON, randomVector.normalized());
                indirectColor = getColor(randomRay, depth - 1) * sphere.albedo;

                // Contribution directe
                for (size_t i = 0; i < lights.size(); i++)
                {
                    Vector axisVector = (P - lights[i].center).normalized();
                    Vector randomVector = (P - lights[i].center).generateRandomCosineVector().normalized();
                    Vector randomLightSource = randomVector * lights[i].radius + lights[i].center;
                    Vector wi = (randomLightSource - P).normalized();
                    double lightDistanceSquared = (randomLightSource - P).norm2();
                    Ray lightRay(P + N * EPSILON, wi);
                    Vector Plight, Nlight;
                    // int objectIndex;
                    // double tlight;
                    // bool intersectLight = this->intersect(lightRay, Plight, Nlight, objectIndex, tlight);

                    // if (intersectLight && tlight * tlight < lightDistanceSquared * 0.99)
                    // {
                    //     continue;
                    // }
                    diffusedColor += lights[i].lightIntensity / (4 * PI * lightDistanceSquared) * std::max(0.0, dot(N, wi)) * dot(randomVector, (-1) * wi) / dot(axisVector, randomVector) * lights[i].albedo * sphere.albedo;
                }
            }

            // réflexion
            if (reflectance > 0)
            {
                Vector reflexionVector = ray.direction - 2 * dot(ray.direction, N) * N;
                Ray reflectedRay(P + N * EPSILON, reflexionVector.normalized());
                Vector reflectedColor = getColor(reflectedRay, depth - 1);

                return (1 - reflectance) * (diffusedColor + indirectColor) + reflectance * reflectedColor;
            }

            return diffusedColor + indirectColor;
        }
    }
    return Vector(0, 0, 0);
}
