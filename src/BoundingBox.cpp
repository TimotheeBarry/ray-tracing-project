#include "../include/BoundingBox.hpp"
#include <limits>
#include <iostream>

bool BoundingBox::intersect(Ray &ray) const
{
    double intervalMin(-std::numeric_limits<double>::infinity()), intervalMax(std::numeric_limits<double>::infinity());
    double tMin, tMax;
    for (int i = 0; i < 3; i++)
    {
        // pour chaque cooordonnées, on calcule l'intervalle d'intersection avec le rayon
        tMin = (min[i] - ray.origin[i]) / ray.direction[i];
        tMax = (max[i] - ray.origin[i]) / ray.direction[i];
        if (tMin > tMax)
        {
            std::swap(tMin, tMax);
        }
        // On met à jour l'intervalle d'intersection global
        if (tMin > intervalMin)
        {
            intervalMin = tMin;
        }
        if (tMax < intervalMax)
        {
            intervalMax = tMax;
        }
        // Si l'intervalle d'intersection est vide, on renvoie false
        if (intervalMin > intervalMax)
        {
            return false;
        }
    }
    return true;
};

void BoundingBox::scale(double s, Vector &center)
{
    min = (min - center) * s + center;
    max = (max - center) * s + center;
};

void BoundingBox::translate(Vector &t)
{
    min = min + t;
    max = max + t;
};

void BoundingBox::computeDimensions(std::vector<Vector> &vertices)
{
    for (int i = 0; i < vertices.size(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (vertices[i][j] < min[j])
            {
                min[j] = vertices[i][j];
            }
            if (vertices[i][j] > max[j])
            {
                max[j] = vertices[i][j];
            }
        }
    }
}
