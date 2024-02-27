#include "../include/Scene.hpp"
#include "../include/Functions.hpp"
#include "../include/Constants.hpp"
#include "../include/Sphere.hpp"
#include "../include/LightSource.hpp"
#include <cmath>
#include <iostream>

void Scene::addObject(const Object &obj)
{
    objects.push_back(&obj);
}

bool Scene::intersect(Ray &ray, Vector &P, Vector &N, int &objectIndex, double &t)
{
    t = 1e99;
    bool intersect = false;
    Vector tempP, tempN;
    for (size_t i = 0; i < objects.size(); i++)
    {
        double tMin = objects[i]->intersect(ray, tempP, tempN);
        if (tMin > 0 && tMin < t)
        {
            t = tMin;
            P = tempP;
            N = tempN;
            objectIndex = i;
            intersect = true;
        }
    }
    return intersect;
}

Vector Scene::getColor(Ray &ray, int depth, bool isIndirect)
{
    if (depth < 0)
    {
        return Vector(0, 0, 0);
    }

    Vector P, N;
    int intersectIndex = -1;
    double t;
    if (this->intersect(ray, P, N, intersectIndex, t))
    {

        const Object *object = objects[intersectIndex];

        // Si on a tapé une source de lumière
        if (LightSource *light = dynamic_cast<LightSource *>(const_cast<Object *>(object)))
        {
            if (isIndirect)
            {
                // Si c'est un rebond de lumière indirecte, on ne prend pas en compte la lumière des sphères lumineuses
                return Vector(0, 0, 0);
            }
            return Vector(1, 1, 1) * light->realIntensity();
        }

        else if (Sphere *sphere = dynamic_cast<Sphere *>(const_cast<Object *>(object)))
        {
            // Si la sphere est transparente
            if (sphere->opacity < 1)
            {
                // entrée ou sortie de la sphère
                bool isEntering = dot(ray.direction, N) < 0;
                double n1, n2;
                if (isEntering)
                {
                    n1 = 1.0;
                    n2 = sphere->refractiveIndex;
                    N = (-1) * N;
                }
                else
                {
                    n1 = sphere->refractiveIndex;
                    n2 = 1.0;
                }
                double k0 = sqr((n1 - n2) / (n1 + n2));

                double R = k0 + (1 - k0) * std::pow(1 - std::abs(dot(ray.direction, N)), 5);
                double T = 1 - R;

                double n = n1 / n2;
                double cosThetaI = dot(ray.direction, N);
                Vector tN = std::sqrt(1 - sqr(n) * (1 - sqr(cosThetaI))) * N;

                Vector tT = n * ray.direction;
                // on choisit aléatoirement entre réflexion et transmission avec proba qui dépend des coefficients R et T
                if (uniform(gen) < T)
                {
                    // transmission
                    Ray ray = Ray(P + N * EPSILON, (tN + tT).normalized());
                    return this->getColor(ray, depth - 1, isIndirect);
                }
                else
                {
                    // réflexion
                    Vector reflexionVector = ray.direction - 2 * dot(ray.direction, N) * N;
                    Ray reflectedRay(P + N * EPSILON, reflexionVector.normalized());
                    return this->getColor(reflectedRay, depth - 1, isIndirect);
                }
            }
            else
            {
                Vector indirectColor(0, 0, 0), diffusedColor(0, 0, 0);
                // on calcule la contribution directe et indirecte que si la sphère n'est pas un miroir pur pour éviter les calculs inutiles
                if (sphere->reflectance < 1)
                {
                    // contribution indirecte
                    Vector randomVector = N.generateRandomCosineVector();
                    Ray randomRay = Ray(P + N * EPSILON, randomVector.normalized());
                    indirectColor = this->getColor(randomRay, depth - 1, true) * sphere->albedo;

                    // Contribution directe
                    for (size_t i = 0; i < objects.size(); i++)
                        if (LightSource *light = dynamic_cast<LightSource *>(const_cast<Object *>(objects[i])))
                        {
                            Vector axisVector = (P - light->center).normalized();
                            Vector randomVector = (P - light->center).generateRandomCosineVector().normalized();
                            Vector randomLightSource = randomVector * light->radius + light->center;
                            Vector wi = (randomLightSource - P).normalized();
                            double lightDistanceSquared = (randomLightSource - P).norm2();
                            Ray lightRay(P + N * EPSILON, wi);
                            Vector Plight, Nlight;
                            int objectIndex;
                            double tlight;
                            if (this->intersect(lightRay, Plight, Nlight, objectIndex, tlight) && tlight * tlight < lightDistanceSquared * (1 - EPSILON))
                            {
                                continue;
                            }
                            diffusedColor += light->realIntensity() / (4 * PI * lightDistanceSquared) * std::max(0.0, dot(N, wi)) * dot(randomVector, (-1) * wi) / dot(axisVector, randomVector) * sphere->albedo;
                        }
                }

                // réflexion
                if (sphere->reflectance > 0)
                {
                    Vector reflexionVector = ray.direction - 2 * dot(ray.direction, N) * N;
                    Ray reflectedRay(P + N * EPSILON, reflexionVector.normalized());
                    Vector reflectedColor = this->getColor(reflectedRay, depth - 1, isIndirect);

                    return (1 - sphere->reflectance) * (diffusedColor + indirectColor) + sphere->reflectance * reflectedColor;
                }

                return diffusedColor + indirectColor;
            }
        }
    }
    return Vector(0, 0, 0);
}
