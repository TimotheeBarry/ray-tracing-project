#include "../include/BVH.hpp"

std::string BVH::printTree(int depth)
{
    std::string s = "";
    for (int i = 0; i < depth; i++)
    {
        s += ".";
    }
    if (isLeaf)
    {
        s += "Leaf: " + std::to_string(iMin) + " " + std::to_string(iMax) + " (" + std::to_string(iMax - iMin + 1) + " triangles)\n";
    }
    else
    {
        s += "Node: " + std::to_string(iMin) + " " + std::to_string(iMax) + "\n";
        if (left != nullptr) {
            s += left->printTree(depth + 1);
        }
        if (right != nullptr) {
            s += right->printTree(depth + 1);
        }
    }
    return s;
}