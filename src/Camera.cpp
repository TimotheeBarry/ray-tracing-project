#include "../include/Camera.hpp"

Ray Camera::launchRay(int i, int j) const
{
    // ouverture de l'objectif
    auto dx = (uniform(gen) - .5) * this->aperture;
    auto dy = (uniform(gen) - .5) * this->aperture;
    Vector rayOrigin = this->origin + Vector(dx, dy, 0);
    // anti aliasing
    double di, dj;
    boxMuller(0.2, di, dj);
    // direction du rayon
    Vector direction = Vector(
        j + .5 - this->width / 2 + dj,
        -(i + .5) + this->height / 2 + di,
        -this->width / (2 * std::tan(this->fov * (PI) / (2 * 180))));
    direction.normalize();
    // point dans le plan focal
    auto target = this->origin + direction * this->focalLength;
    return Ray(rayOrigin, (target - rayOrigin).normalized());
}