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
            return computeColorFromBRDF(ray, P, N, depth, isIndirect, albedo, sphere->mat);
        }
        else if (TriangleMesh *mesh = dynamic_cast<TriangleMesh *>(const_cast<Object *>(object)))
        {
            return computeColorFromBRDF(ray, P, N, depth, isIndirect, albedo, mesh->mat);
        }
    }
    return Vector(0, 0, 0);
}

Vector Scene::getIndirectColor(Ray &ray, Vector &P, Vector &N, int depth, Vector &albedo, Material *mat)
{
    Vector randomVector;
    if (BlinnPhong *material = dynamic_cast<BlinnPhong *>(const_cast<Material *>(mat)))
    {
        double p = 1 - material->shininess;
        if (uniform(gen) < p)
        {
            // Diffusion
            randomVector = generateRandomCosineVector(N);
        }
        else
        {
            // Réflexion spéculaire (Blinn-Phong)
            Vector reflexionVector = ray.direction - 2 * dot(ray.direction, N) * N;
            randomVector = generateRandomBlinnPhongVector(reflexionVector, material->alpha);
            if (dot(randomVector, N) < 0)
            {
                // Si le vecteur aléatoire va sous la surface, on le rejette
                return Vector(0, 0, 0);
            }
        }
        double proba = p * dot(randomVector, N) / PI + (1 - p) * (material->alpha + 1) * std::pow(dot(randomVector, ray.direction), material->alpha) / (2 * PI);
        Ray randomRay = Ray(P + N * EPSILON, randomVector.normalized());
        return this->getColor(randomRay, depth, true) * material->brdf(albedo, randomVector, ray.direction, N) * dot(randomVector, N) / proba;
    }
    else
    {
        randomVector = generateRandomCosineVector(N);
        Ray randomRay = Ray(P + N * EPSILON, randomVector.normalized());
        return this->getColor(randomRay, depth, true) * albedo;
    }
}

Vector Scene::getTransparentColor(Ray &ray, Vector &P, Vector &N, int depth, bool isIndirect, double refractiveIndex)
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
        return getReflectionColor(ray, P, N, depth - 1, isIndirect);
    }
}

Vector Scene::getDiffusedColor(Ray &ray, Vector &P, Vector &N, int depth, bool isIndirect, Vector &albedo, Material *mat)
{
    Vector color(0, 0, 0);
    for (size_t i = 0; i < objects.size(); i++)
    {
        if (LightSource *light = dynamic_cast<LightSource *>(const_cast<Object *>(objects[i])))
        {
            Vector axisVector = (P - light->center);
            Vector randomVector = generateRandomCosineVector(axisVector).normalized();
            axisVector.normalize();
            Vector randomLightSource = randomVector * light->radius + light->center;
            Vector wi = (randomLightSource - P).normalized();
            double lightDistanceSquared = (randomLightSource - P).normSquared();
            Ray lightRay(P + N * EPSILON, wi, sqrt(lightDistanceSquared));
            if (this->intersectObjectOnly(lightRay))
            {
                continue;
            }
            color += light->realIntensity() / (4 * PI * lightDistanceSquared) * std::max(0.0, dot(N, wi)) * dot(randomVector, (-1) * wi) / dot(axisVector, randomVector) * mat->brdf(albedo, wi, ray.direction, N);
        }
    }
    color += getIndirectColor(ray, P, N, depth - 1, albedo, mat);
    return color;
}

Vector Scene::getReflectionColor(Ray &ray, Vector &P, Vector &N, int depth, bool isIndirect)
{
    Vector reflexionVector = ray.direction - 2 * dot(ray.direction, N) * N;
    Ray reflectedRay(P + N * EPSILON, reflexionVector.normalized());
    return this->getColor(reflectedRay, depth, isIndirect);
}

Vector Scene::getMirrorColor(Ray &ray, Vector &P, Vector &N, int depth, bool isIndirect, double reflectance, Vector &albedo, Material *mat)
{
    Vector diffused(0, 0, 0);
    if (reflectance < 1)
    {
        // pas la peine de calculer la couleur diffusée si la sphère est totalement réfléchissante
        diffused = getDiffusedColor(ray, P, N, depth, isIndirect, albedo, mat);
    }
    Vector reflected = getReflectionColor(ray, P, N, depth - 1, isIndirect);

    return (1 - reflectance) * (diffused) + reflectance * reflected;
}

Vector Scene::computeColorFromBRDF(Ray &ray, Vector &P, Vector &N, int depth, bool isIndirect, Vector &albedo, Material *mat)
{
    if (Transparent *material = dynamic_cast<Transparent *>(const_cast<Material *>(mat)))
    {
        return getTransparentColor(ray, P, N, depth, isIndirect, material->refractiveIndex);
    }
    if (Mirror *material = dynamic_cast<Mirror *>(const_cast<Material *>(mat)))
    {
        return getMirrorColor(ray, P, N, depth, isIndirect, material->reflectance, albedo, mat);
    }
    else if (Diffuse *material = dynamic_cast<Diffuse *>(const_cast<Material *>(mat)))
    {
        return getDiffusedColor(ray, P, N, depth, isIndirect, albedo, mat);
    }
    else if (BlinnPhong *material = dynamic_cast<BlinnPhong *>(const_cast<Material *>(mat)))
    {
        return getDiffusedColor(ray, P, N, depth, isIndirect, albedo, mat);
    }
    return Vector(0, 0, 0);
}
