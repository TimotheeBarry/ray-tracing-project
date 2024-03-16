#include "../include/Scene.hpp"
#include "../include/Functions.hpp"
#include "../include/Constants.hpp"
#include "../include/Sphere.hpp"
#include "../include/LightSource.hpp"
#include "../include/TriangleMesh.hpp"
#include <cmath>
#include <iostream>

void Scene::addObject(const Object &obj)
{
    objects.push_back(&obj);
}

bool Scene::intersect(Ray &ray, Vector &P, Vector &N, int &objectIndex, double &t, Vector &albedo)
{
    t = 1e99;
    bool intersect = false;
    Vector tempP, tempN, tempAlbedo;
    for (size_t i = 0; i < objects.size(); i++)
    {
        double tMin = objects[i]->intersect(ray, tempP, tempN, tempAlbedo);
        if (tMin > 0 && tMin < t)
        {
            t = tMin;
            P = tempP;
            N = tempN;
            albedo = tempAlbedo;
            objectIndex = i;
            intersect = true;
        }
    }
    return intersect;
}

// méthode qui renvoie true si le rayon intersecte un objet, sans calculer le point d'intersection
// utilisée pour les ombres (plus rapide)
bool Scene::intersectObjectOnly(Ray &ray)
{
    for (size_t i = 0; i < objects.size(); i++)
    {
        // si l'objet est une lightsource, on ne la compte pas
        if (LightSource *light = dynamic_cast<LightSource *>(const_cast<Object *>(objects[i])))
        {
            continue;
        }
        if (objects[i]->fastIntersect(ray))
        {
            return true;
        }
    }
    return false;
}

Vector Scene::getColor(Ray &ray, int depth, bool isIndirect)
{
    if (depth < 0)
    {
        return Vector(0, 0, 0);
    }

    Vector P, N, albedo;
    int objectIndex = -1;
    double t;
    if (this->intersect(ray, P, N, objectIndex, t, albedo))
    {
        const Object *object = objects[objectIndex];

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
            if (Transparent *brdf = dynamic_cast<Transparent *>(const_cast<BRDF *>(sphere->brdf)))
            {
                return getTransparentSphereColor(ray, P, N, depth, isIndirect, brdf->refractiveIndex);
            }
            else
            {
                if (Mirror *brdf = dynamic_cast<Mirror *>(const_cast<BRDF *>(sphere->brdf)))
                {
                    // contribution from reflexion
                    Vector reflexion = getReflectionColor(ray, P, N, depth, isIndirect);
                    // contribution from diffusion (direct and indirect)
                    Vector diffused = getDiffusedColor(ray, P, N, depth, isIndirect, albedo);
                    return (1 - brdf->reflectance) * (diffused + reflexion) + brdf->reflectance * reflexion;
                }
                else if (LambertianBRDF *brdf = dynamic_cast<LambertianBRDF *>(const_cast<BRDF *>(sphere->brdf)))
                {
                    return getDiffusedColor(ray, P, N, depth, isIndirect, albedo);
                }
                else if (BlinnPhongBRDF *brdf = dynamic_cast<BlinnPhongBRDF *>(const_cast<BRDF *>(sphere->brdf)))
                {
                    // TODO
                }
            }
        }
        else if (TriangleMesh *mesh = dynamic_cast<TriangleMesh *>(const_cast<Object *>(object)))
        {

            Vector color = getDiffusedColor(ray, P, N, depth, isIndirect, albedo);
            // réflexion
            if (mesh->reflectance > 0)
            {
                Vector reflectedColor = getReflectionColor(ray, P, N, depth, isIndirect);

                return (1 - mesh->reflectance) * color + mesh->reflectance * reflectedColor;
            }
            return color;
        }
    }
    return Vector(0, 0, 0);
}

Vector Scene::getIndirectColor(Vector &P, Vector &N, int depth, Vector &albedo)
{
    Vector randomVector = N.generateRandomCosineVector();
    Ray randomRay = Ray(P + N * EPSILON, randomVector.normalized());
    return this->getColor(randomRay, depth - 1, true) * albedo;
}

Vector Scene::getTransparentSphereColor(Ray &ray, Vector &P, Vector &N, int depth, bool isIndirect, double refractiveIndex)
{
    // entrée ou sortie de la sphère
    bool isEntering = dot(ray.direction, N) < 0;
    double n1, n2;
    if (isEntering)
    {
        n1 = 1.0;
        n2 = refractiveIndex;
        N = (-1) * N;
    }
    else
    {
        n1 = refractiveIndex;
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

Vector Scene::getDiffusedColor(Ray &ray, Vector &P, Vector &N, int depth, bool isIndirect, Vector &albedo)
{
    Vector color(0, 0, 0);
    for (size_t i = 0; i < objects.size(); i++)
    {
        if (LightSource *light = dynamic_cast<LightSource *>(const_cast<Object *>(objects[i])))
        {
            Vector axisVector = (P - light->center).normalized();
            Vector randomVector = (P - light->center).generateRandomCosineVector().normalized();
            Vector randomLightSource = randomVector * light->radius + light->center;
            Vector wi = (randomLightSource - P).normalized();
            double lightDistanceSquared = (randomLightSource - P).normSquared();
            Ray lightRay(P + N * EPSILON, wi, sqrt(lightDistanceSquared));
            if (this->intersectObjectOnly(lightRay))
            {
                continue;
            }
            color += light->realIntensity() / (4 * PI * lightDistanceSquared) * std::max(0.0, dot(N, wi)) * dot(randomVector, (-1) * wi) / dot(axisVector, randomVector) * albedo;
        }
    }
    return color + getIndirectColor(P, N, depth, albedo);
}

Vector Scene::getReflectionColor(Ray &ray, Vector &P, Vector &N, int depth, bool isIndirect)
{
    Vector reflexionVector = ray.direction - 2 * dot(ray.direction, N) * N;
    Ray reflectedRay(P + N * EPSILON, reflexionVector.normalized());
    return this->getColor(reflectedRay, depth - 1, isIndirect);
}
