#pragma once

#include <vector>
#include "Object.hpp"
#include "BoundingBox.hpp"
#include "TriangleIndices.hpp"
#include "BVH.hpp"

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
    BVH bvh;
    int bvhMinTriangles = 5;

    void readOBJ(const char *obj);
    void scale(double s);
    void translate(Vector t);
    void rotate(double angle, Vector axis);
    double intersect(Ray &ray, Vector &P, Vector &N) const override;
    bool fastIntersect(Ray &ray) const override;
    Vector getBarycenter() const;
    void updateBoundingBox(BoundingBox &bbox, int iMin, int iMax);
    void updateMainBoundingBox();
    double getTriangleCenterAlongAxis(int index, int axis) const;
    void buildBVH(BVH &node, int iMin, int iMax);
    void initBVH();
};