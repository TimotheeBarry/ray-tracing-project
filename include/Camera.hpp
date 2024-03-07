#include "Vector.hpp"
#include "Ray.hpp"
#include "Constants.hpp"
#include "Functions.hpp"

class Camera
{
public:
    Camera(
        Vector origin,
        int width,
        int height,
        double fov,
        double aperture,
        double focalLength) : origin(origin),
                              fov(fov), width(width),
                              height(height),
                              aperture(aperture),
                              focalLength(focalLength) {}

    Vector origin;
    int width;
    int height;
    double fov;
    double aperture;
    double focalLength;

    Ray launchRay(int i, int j) const;
};