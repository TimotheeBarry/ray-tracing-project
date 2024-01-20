#pragma once
#include "Vector.hpp"

class Ray
{
public:
    // explicit Ray(Vector O, Vector u);
    Ray(Vector origin, Vector direction);

    Vector origin;
    Vector direction;
};