#pragma once

#include <vector>
#include "Object.hpp"
#include "BoundingBox.hpp"

class TriangleIndices
{
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group){};
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;    // indices within the uv coordinates array
    int ni, nj, nk;       // indices within the normals array
    int group;            // face group
};

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

    void readOBJ(const char *obj);
    void scale(double s);
    void translate(Vector t);
    void rotate(double angle, Vector axis);
    double intersect(Ray &ray, Vector &P, Vector &N) const override;
    bool fastIntersect(Ray &ray) const override;
    Vector getBarycenter() const;
    BoundingBox getBoundingBox() const;
    BoundingBox bbox;

private:
        void computeBoundingBox();
};