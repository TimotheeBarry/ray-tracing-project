#pragma once

#include <vector>
#include "Object.hpp"
#include "BoundingBox.hpp"
#include "TriangleIndices.hpp"

class TriangleMesh : public Object
{
public:
    ~TriangleMesh() {}
    TriangleMesh(){};

    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    BoundingBox bbox;
    int bboxMaxDepth = 3;
    int bboxMinTriangles = 2;

    void readOBJ(const char *obj);
    void scale(double s);
    void translate(Vector t);
    void rotate(double angle, Vector axis);
    double intersect(Ray &ray, Vector &P, Vector &N) const override;
    bool fastIntersect(Ray &ray) const override;
    Vector getBarycenter() const;
    void updateBoundingBox();
    void initializeBoundingBoxDimensions();
};