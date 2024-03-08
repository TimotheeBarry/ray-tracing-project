#pragma once

#include "BoundingBox.hpp"

class BVH
{
public:
    BVH() : left(nullptr), right(nullptr), iMin(0), iMax(0), isLeaf(false) {}
    ~BVH()
    {
        delete left;
        delete right;
    }

    BVH *left;
    BVH *right;
    int iMin, iMax;
    bool isLeaf;
    BoundingBox bbox;

    std::string printTree(int depth = 0);
};