#pragma once

#include <vector>
#include "Object.hpp"
#include "BoundingBox.hpp"
#include "TriangleIndices.hpp"
#include "BVH.hpp"
#include "Material.hpp"
#include "Constants.hpp"

class TriangleMesh : public Object
{
public:
    ~TriangleMesh() {}
    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    BoundingBox bbox;
    BVH bvh;
    Material *mat;
    int bvhMinTriangles;
    
    TriangleMesh(const char *obj,
                 const char *textures,
                 Material *mat = nullptr,
                 int bvhMinTriangles = 5)
        : mat(mat), bvhMinTriangles(bvhMinTriangles)
    {
        readOBJ(obj);
        readPNGTexture(textures);
        if (mat == nullptr)
        {
            this->mat = new Diffuse();
        }
    }

        void readOBJ(const char *obj);
    void readPNGTexture(const char *filename);
    void scale(double s);
    void translate(Vector t);
    void rotate(double angle, Vector axis);
    double intersect(Ray &ray, Vector &P, Vector &N, Vector &albedo) const override;
    bool fastIntersect(Ray &ray) const override;
    Vector getBarycenter() const;
    void initBVH();

private:
    std::vector<std::vector<unsigned char>> textures;
    std::vector<int> w, h;
    double getTriangleCenterAlongAxis(int index, int axis) const;
    void updateBoundingBox(BoundingBox &bbox, int iMin, int iMax);
    void buildBVH(BVH &node, int iMin, int iMax);
    void updateMainBoundingBox();
};