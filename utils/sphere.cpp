#include "vector.cpp"

class Sphere {
public:
    explicit Sphere(Vector C, double R)
	{
		center = C;
        radius = R;
    }

    Vector center;
    double radius;
};