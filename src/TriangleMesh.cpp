#include "../include/TriangleMesh.hpp"
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <stack>
#include <cmath>
#include "../stb/stb_image.h"

double TriangleMesh::intersect(Ray &ray, Vector &P, Vector &N, Vector &albedo) const
{
    // vérifie premièrement si le rayon intersecte la bounding box
    if (!bvh.bbox.intersect(ray))
    {
        return -1;
    }
    else
    {
        double tmin = 1e99;
        // parcours en profondeur de l'arbre
        std::stack<const BVH *> stack;
        stack.push(&bvh);
        while (!stack.empty())
        {
            const BoundingBox &bbox = stack.top()->bbox;
            const BVH *node = stack.top();
            stack.pop();
            if (node->isLeaf)
            {
                for (int i = node->iMin; i <= node->iMax; i++)
                {
                    const TriangleIndices &triangle = indices[i];
                    // sommets du triangle
                    Vector v0 = vertices[triangle.vtxi];
                    Vector v1 = vertices[triangle.vtxj];
                    Vector v2 = vertices[triangle.vtxk];

                    Vector e1 = v1 - v0;
                    Vector e2 = v2 - v0;
                    // normale
                    Vector Ni = cross(e1, e2);
                    // calcul de l'intersection
                    double divisor = dot(ray.direction, Ni);
                    if (divisor == 0)
                    {
                        // rayon parallèle au triangle
                        continue;
                    }
                    double t = dot(v0 - ray.origin, Ni) / divisor;

                    if (t < 0 || t > ray.length)
                    {
                        // point d'intersection derrière le rayon
                        continue;
                    }
                    Vector crossNumerator = cross(v0 - ray.origin, ray.direction);
                    double beta = dot(e2, crossNumerator) / divisor;
                    double gamma = -dot(e1, crossNumerator) / divisor;

                    // conditions d'intersection
                    if (beta >= 0 && gamma >= 0 && beta + gamma <= 1 && t < tmin)
                    {
                        P = ray.origin + ray.direction * t;
                        // lissage de phong
                        double alpha = 1 - beta - gamma;
                        if (triangle.ni == -1 || triangle.nj == -1 || triangle.nk == -1)
                        {
                            N = Ni;
                        }
                        else
                        {
                            N = normals[triangle.ni] * alpha + normals[triangle.nj] * beta + normals[triangle.nk] * gamma;
                        }
                        tmin = t;
                        // texture
                        if (!textures.empty() && !textures[triangle.group].empty() && triangle.uvi != -1 && triangle.uvj != -1 && triangle.uvk != -1)
                        {
                            Vector uv = uvs[triangle.uvi] * alpha + uvs[triangle.uvj] * beta + uvs[triangle.uvk] * gamma;

                            // récupération de la couleur du pixel (interpolation bilinéaire
                            int w = this->w[triangle.group];
                            int h = this->h[triangle.group];
                            int x = uv[0] * (w - 1);
                            int y = h - 1 - uv[1] * h;
                            int index = 3 * (y * w + x);
                            double r = textures[triangle.group][index] / 255.0;
                            double g = textures[triangle.group][index + 1] / 255.0;
                            double b = textures[triangle.group][index + 2] / 255.0;
                            albedo = pow(Vector(r, g, b), GAMMA);
                        }
                        else if (!vertexcolors.empty())
                        {
                            albedo = vertexcolors[triangle.vtxi] * alpha + vertexcolors[triangle.vtxj] * beta + vertexcolors[triangle.vtxk] * gamma;
                        }
                        else
                        {
                            albedo = Vector(1, 1, 1);
                        }
                    }
                }
            }
            else
            {
                double t1, t2;
                if (node->left->bbox.intersect(ray))
                {
                    stack.push(node->left);
                }
                if (node->right->bbox.intersect(ray))
                {
                    stack.push(node->right);
                }
            }
        }
        return tmin;
    }
}

bool TriangleMesh::fastIntersect(Ray &ray) const
{
    // vérifie premièrement si le rayon intersecte la bounding box
    if (!bvh.bbox.intersect(ray))
    {
        return false;
    }
    else
    {
        // parcours en profondeur de l'arbre
        std::stack<const BVH *> stack;
        stack.push(&bvh);
        while (!stack.empty())
        {
            const BoundingBox &bbox = stack.top()->bbox;
            const BVH *node = stack.top();
            stack.pop();
            if (node->isLeaf)
            {
                for (int i = node->iMin; i <= node->iMax; i++)
                {
                    const TriangleIndices &triangle = indices[i];
                    // sommets du triangle
                    Vector v0 = vertices[triangle.vtxi];
                    Vector v1 = vertices[triangle.vtxj];
                    Vector v2 = vertices[triangle.vtxk];

                    Vector e1 = v1 - v0;
                    Vector e2 = v2 - v0;
                    // normale
                    Vector Ni = cross(e1, e2);
                    // calcul de l'intersection
                    double divisor = dot(ray.direction, Ni);
                    if (divisor == 0)
                    {
                        // rayon parallèle au triangle
                        continue;
                    }
                    double t = dot(v0 - ray.origin, Ni) / divisor;

                    if (t < 0 || t > ray.length)
                    {
                        // point d'intersection derrière le rayon
                        continue;
                    }
                    Vector crossNumerator = cross(v0 - ray.origin, ray.direction);
                    double beta = dot(e2, crossNumerator) / divisor;
                    double gamma = -dot(e1, crossNumerator) / divisor;

                    // conditions d'intersection
                    if (beta >= 0 && gamma >= 0 && beta + gamma <= 1)
                    {
                        return true;
                    }
                }
            }
            else
            {
                double t1, t2;
                if (node->left->bbox.intersect(ray))
                {
                    stack.push(node->left);
                }
                if (node->right->bbox.intersect(ray))
                {
                    stack.push(node->right);
                }
            }
        }
        return false;
    }
}

Vector TriangleMesh::getBarycenter() const
{
    Vector barycenter;
    for (int i = 0; i < vertices.size(); i++)
    {
        barycenter = barycenter + vertices[i];
    }
    barycenter = barycenter / vertices.size();
    return barycenter;
}

void TriangleMesh::translate(Vector t)
{
    for (int i = 0; i < vertices.size(); i++)
    {
        vertices[i] = vertices[i] + t;
    }
    bbox.translate(t);
}

void TriangleMesh::rotate(double angle, Vector axis)
{
    Vector barycenter = getBarycenter();
    for (int i = 0; i < vertices.size(); i++)
    {
        vertices[i] = barycenter + (vertices[i] - barycenter).rotate(angle, axis);
    }
    updateMainBoundingBox();
}

void TriangleMesh::scale(double s)
{
    Vector barycenter = getBarycenter();
    for (int i = 0; i < vertices.size(); i++)
    {
        vertices[i] = barycenter + (vertices[i] - barycenter) * s;
    }
    bbox.scale(s, barycenter);
}

void TriangleMesh::readOBJ(const char *obj)
{

    char matfile[255];
    char grp[255];

    FILE *f;
    f = fopen(obj, "r");
    int curGroup = -1;
    while (!feof(f))
    {
        char line[255];
        if (!fgets(line, 255, f))
            break;

        std::string linetrim(line);
        linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
        strcpy(line, linetrim.c_str());

        if (line[0] == 'u' && line[1] == 's')
        {
            sscanf(line, "usemtl %[^\n]\n", grp);
            curGroup++;
        }

        if (line[0] == 'v' && line[1] == ' ')
        {
            Vector vec;

            Vector col;
            if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6)
            {
                col[0] = std::min(1., std::max(0., col[0]));
                col[1] = std::min(1., std::max(0., col[1]));
                col[2] = std::min(1., std::max(0., col[2]));

                vertices.push_back(vec);
                vertexcolors.push_back(col);
            }
            else
            {
                sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                vertices.push_back(vec);
            }
        }
        if (line[0] == 'v' && line[1] == 'n')
        {
            Vector vec;
            sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
            normals.push_back(vec);
        }
        if (line[0] == 'v' && line[1] == 't')
        {
            // Vector vec;
            // sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
            // uvs.push_back(vec);

            Vector vec;
            double u, v;
            sscanf(line, "vt %lf %lf\n", &u, &v);
            // Ensure UV coordinates fall within the range [0, 1]
            u -= floor(u);
            v -= floor(v);
            vec[0] = u;
            vec[1] = v;
            uvs.push_back(vec);
        }
        if (line[0] == 'f')
        {
            TriangleIndices t;
            int i0, i1, i2, i3;
            int j0, j1, j2, j3;
            int k0, k1, k2, k3;
            int nn;
            t.group = curGroup;

            char *consumedline = line + 1;
            int offset;

            nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
            if (nn == 9)
            {
                if (i0 < 0)
                    t.vtxi = vertices.size() + i0;
                else
                    t.vtxi = i0 - 1;
                if (i1 < 0)
                    t.vtxj = vertices.size() + i1;
                else
                    t.vtxj = i1 - 1;
                if (i2 < 0)
                    t.vtxk = vertices.size() + i2;
                else
                    t.vtxk = i2 - 1;
                if (j0 < 0)
                    t.uvi = uvs.size() + j0;
                else
                    t.uvi = j0 - 1;
                if (j1 < 0)
                    t.uvj = uvs.size() + j1;
                else
                    t.uvj = j1 - 1;
                if (j2 < 0)
                    t.uvk = uvs.size() + j2;
                else
                    t.uvk = j2 - 1;
                if (k0 < 0)
                    t.ni = normals.size() + k0;
                else
                    t.ni = k0 - 1;
                if (k1 < 0)
                    t.nj = normals.size() + k1;
                else
                    t.nj = k1 - 1;
                if (k2 < 0)
                    t.nk = normals.size() + k2;
                else
                    t.nk = k2 - 1;
                indices.push_back(t);
            }
            else
            {
                nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                if (nn == 6)
                {
                    if (i0 < 0)
                        t.vtxi = vertices.size() + i0;
                    else
                        t.vtxi = i0 - 1;
                    if (i1 < 0)
                        t.vtxj = vertices.size() + i1;
                    else
                        t.vtxj = i1 - 1;
                    if (i2 < 0)
                        t.vtxk = vertices.size() + i2;
                    else
                        t.vtxk = i2 - 1;
                    if (j0 < 0)
                        t.uvi = uvs.size() + j0;
                    else
                        t.uvi = j0 - 1;
                    if (j1 < 0)
                        t.uvj = uvs.size() + j1;
                    else
                        t.uvj = j1 - 1;
                    if (j2 < 0)
                        t.uvk = uvs.size() + j2;
                    else
                        t.uvk = j2 - 1;
                    indices.push_back(t);
                }
                else
                {
                    nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                    if (nn == 3)
                    {
                        if (i0 < 0)
                            t.vtxi = vertices.size() + i0;
                        else
                            t.vtxi = i0 - 1;
                        if (i1 < 0)
                            t.vtxj = vertices.size() + i1;
                        else
                            t.vtxj = i1 - 1;
                        if (i2 < 0)
                            t.vtxk = vertices.size() + i2;
                        else
                            t.vtxk = i2 - 1;
                        indices.push_back(t);
                    }
                    else
                    {
                        nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                        if (i0 < 0)
                            t.vtxi = vertices.size() + i0;
                        else
                            t.vtxi = i0 - 1;
                        if (i1 < 0)
                            t.vtxj = vertices.size() + i1;
                        else
                            t.vtxj = i1 - 1;
                        if (i2 < 0)
                            t.vtxk = vertices.size() + i2;
                        else
                            t.vtxk = i2 - 1;
                        if (k0 < 0)
                            t.ni = normals.size() + k0;
                        else
                            t.ni = k0 - 1;
                        if (k1 < 0)
                            t.nj = normals.size() + k1;
                        else
                            t.nj = k1 - 1;
                        if (k2 < 0)
                            t.nk = normals.size() + k2;
                        else
                            t.nk = k2 - 1;
                        indices.push_back(t);
                    }
                }
            }

            consumedline = consumedline + offset;

            while (true)
            {
                if (consumedline[0] == '\n')
                    break;
                if (consumedline[0] == '\0')
                    break;
                nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                TriangleIndices t2;
                t2.group = curGroup;
                if (nn == 3)
                {
                    if (i0 < 0)
                        t2.vtxi = vertices.size() + i0;
                    else
                        t2.vtxi = i0 - 1;
                    if (i2 < 0)
                        t2.vtxj = vertices.size() + i2;
                    else
                        t2.vtxj = i2 - 1;
                    if (i3 < 0)
                        t2.vtxk = vertices.size() + i3;
                    else
                        t2.vtxk = i3 - 1;
                    if (j0 < 0)
                        t2.uvi = uvs.size() + j0;
                    else
                        t2.uvi = j0 - 1;
                    if (j2 < 0)
                        t2.uvj = uvs.size() + j2;
                    else
                        t2.uvj = j2 - 1;
                    if (j3 < 0)
                        t2.uvk = uvs.size() + j3;
                    else
                        t2.uvk = j3 - 1;
                    if (k0 < 0)
                        t2.ni = normals.size() + k0;
                    else
                        t2.ni = k0 - 1;
                    if (k2 < 0)
                        t2.nj = normals.size() + k2;
                    else
                        t2.nj = k2 - 1;
                    if (k3 < 0)
                        t2.nk = normals.size() + k3;
                    else
                        t2.nk = k3 - 1;
                    indices.push_back(t2);
                    consumedline = consumedline + offset;
                    i2 = i3;
                    j2 = j3;
                    k2 = k3;
                }
                else
                {
                    nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                    if (nn == 2)
                    {
                        if (i0 < 0)
                            t2.vtxi = vertices.size() + i0;
                        else
                            t2.vtxi = i0 - 1;
                        if (i2 < 0)
                            t2.vtxj = vertices.size() + i2;
                        else
                            t2.vtxj = i2 - 1;
                        if (i3 < 0)
                            t2.vtxk = vertices.size() + i3;
                        else
                            t2.vtxk = i3 - 1;
                        if (j0 < 0)
                            t2.uvi = uvs.size() + j0;
                        else
                            t2.uvi = j0 - 1;
                        if (j2 < 0)
                            t2.uvj = uvs.size() + j2;
                        else
                            t2.uvj = j2 - 1;
                        if (j3 < 0)
                            t2.uvk = uvs.size() + j3;
                        else
                            t2.uvk = j3 - 1;
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        indices.push_back(t2);
                    }
                    else
                    {
                        nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                        if (nn == 2)
                        {
                            if (i0 < 0)
                                t2.vtxi = vertices.size() + i0;
                            else
                                t2.vtxi = i0 - 1;
                            if (i2 < 0)
                                t2.vtxj = vertices.size() + i2;
                            else
                                t2.vtxj = i2 - 1;
                            if (i3 < 0)
                                t2.vtxk = vertices.size() + i3;
                            else
                                t2.vtxk = i3 - 1;
                            if (k0 < 0)
                                t2.ni = normals.size() + k0;
                            else
                                t2.ni = k0 - 1;
                            if (k2 < 0)
                                t2.nj = normals.size() + k2;
                            else
                                t2.nj = k2 - 1;
                            if (k3 < 0)
                                t2.nk = normals.size() + k3;
                            else
                                t2.nk = k3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            k2 = k3;
                            indices.push_back(t2);
                        }
                        else
                        {
                            nn = sscanf(consumedline, "%u%n", &i3, &offset);
                            if (nn == 1)
                            {
                                if (i0 < 0)
                                    t2.vtxi = vertices.size() + i0;
                                else
                                    t2.vtxi = i0 - 1;
                                if (i2 < 0)
                                    t2.vtxj = vertices.size() + i2;
                                else
                                    t2.vtxj = i2 - 1;
                                if (i3 < 0)
                                    t2.vtxk = vertices.size() + i3;
                                else
                                    t2.vtxk = i3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                indices.push_back(t2);
                            }
                            else
                            {
                                consumedline = consumedline + 1;
                            }
                        }
                    }
                }
            }
        }
    }
    fclose(f);
    updateMainBoundingBox();
};

void TriangleMesh::readPNGTexture(const char *filename)
{
    // int width, height, numChannels;
    // unsigned char *data = stbi_load(filename, &width, &height, &numChannels, 3);

    // if (!data || data == NULL)
    // {
    //     std::cerr << "Error: Couldn't open " << filename << " for reading.\n";
    //     return;
    // }

    // std::vector<unsigned char> texture;

    // // if 4 channels, we have to remove the alpha channel
    // if (numChannels == 4)
    // {
    //     std::vector<unsigned char> texture;
    //     int size = width * height * 3;
    //     texture.resize(size);
    //     for (int i = 0; i < width * height; i += 1)
    //     {
    //         texture[i * 3] = data[i * 4];
    //         texture[i * 3 + 1] = data[i * 4 + 1];
    //         texture[i * 3 + 2] = data[i * 4 + 2];
    //     }
    //     textures.push_back(texture);
    // }
    // else
    // {
    //     std::vector<unsigned char> texture;
    //     int size = width * height * 3;
    //     texture.resize(size);
    //     for (int i = 0; i < width * height * numChannels; i += numChannels)
    //     {
    //         texture[i] = data[i];
    //         texture[i + 1] = data[i + 1];
    //         texture[i + 2] = data[i + 2];
    //     }
    //     textures.push_back(texture);
    // }
    // stbi_image_free(data);

    // textures.push_back(texture);
    // w.push_back(width);
    // h.push_back(height);

    int channels, width, height;

    unsigned char *data = stbi_load(filename, &width, &height, &channels, STBI_rgb);
    if (!data)
    {
        // Handle error if image loading fails
        throw std::runtime_error("Failed to load image.");
    }

    size_t dataSize = width * height * 3;

    std::vector<unsigned char> texture(data, data + dataSize);
    stbi_image_free(data);

    textures.push_back(texture);
    w.push_back(width);
    h.push_back(height);
}

void TriangleMesh::updateBoundingBox(BoundingBox &bbox, int iMin, int iMax)
{
    bbox.min = Vector(1e99, 1e99, 1e99);
    bbox.max = Vector(-1e99, -1e99, -1e99);
    for (int i = iMin; i < iMax; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (vertices[indices[i].vtxi][j] < bbox.min[j])
            {
                bbox.min[j] = vertices[indices[i].vtxi][j];
            }
            if (vertices[indices[i].vtxi][j] > bbox.max[j])
            {
                bbox.max[j] = vertices[indices[i].vtxi][j];
            }
            if (vertices[indices[i].vtxj][j] < bbox.min[j])
            {
                bbox.min[j] = vertices[indices[i].vtxj][j];
            }
            if (vertices[indices[i].vtxj][j] > bbox.max[j])
            {
                bbox.max[j] = vertices[indices[i].vtxj][j];
            }
            if (vertices[indices[i].vtxk][j] < bbox.min[j])
            {
                bbox.min[j] = vertices[indices[i].vtxk][j];
            }
            if (vertices[indices[i].vtxk][j] > bbox.max[j])
            {
                bbox.max[j] = vertices[indices[i].vtxk][j];
            }
        }
    }
}

void TriangleMesh::updateMainBoundingBox()
{
    for (int i = 0; i < vertices.size(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (vertices[i][j] < bbox.min[j])
            {
                bbox.min[j] = vertices[i][j];
            }
            if (vertices[i][j] > bbox.max[j])
            {
                bbox.max[j] = vertices[i][j];
            }
        }
    }
}

double TriangleMesh::getTriangleCenterAlongAxis(int index, int axis) const
{
    const TriangleIndices &triangle = indices[index];
    return (vertices[triangle.vtxi][axis] + vertices[triangle.vtxj][axis] + vertices[triangle.vtxk][axis]) / 3;
}

// construit la BVH avec start inclus, end exclus (indices des triangles)
void TriangleMesh::buildBVH(BVH &node, int start, int end)
{
    // on met à jour le noeud
    node.iMin = start;
    node.iMax = end - 1;
    updateBoundingBox(node.bbox, start, end);

    // si on a le nombre minimal de triangles par bvh, on arrête
    if (end - start <= bvhMinTriangles)
    {
        node.isLeaf = true;
        return;
    }

    // Recherche de la dimension la plus étendue de la bbox
    Vector diag = node.bbox.max - node.bbox.min;
    int dim = 0;
    for (int i = 1; i < 3; i++)
    {
        if (diag[i] > diag[dim])
        {
            dim = i;
        }
    }
    // Définition du pivot
    double leftPointer = node.iMin;
    double rightPointer = node.iMax - 1;
    double leftCenter = getTriangleCenterAlongAxis(leftPointer, dim);
    double rightCenter = getTriangleCenterAlongAxis(rightPointer, dim);
    // algorithm to sort the triangles along the axis with the pivot
    while (leftPointer < rightPointer)
    {
        // Si les deux triangles sont bien placé, on avance le pointeur gauche
        while (leftPointer < rightPointer && leftCenter <= rightCenter)
        {
            leftPointer++;
            leftCenter = getTriangleCenterAlongAxis(leftPointer, dim);
        }
        // Si les deux triangles mal placés, on recule le pointeur droit et on swap
        while (leftPointer < rightPointer && rightCenter < leftCenter)
        {
            std::swap(indices[leftPointer], indices[rightPointer]);
            rightPointer--;
            rightCenter = getTriangleCenterAlongAxis(rightPointer, dim);
        }
    }
    node.left = new BVH();
    node.right = new BVH();

    buildBVH(*node.left, start, leftPointer);
    buildBVH(*node.right, leftPointer, end);
}

void TriangleMesh::initBVH()
{
    bvh = BVH();
    buildBVH(bvh, 0, indices.size());
}