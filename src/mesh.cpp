//
// Created by 胡 on 2023/6/21.
//

#include "../include/mesh.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <utility>
#include <sstream>
#include "../include/box.h"

double Mesh::intersect(const Ray &r, double tmin, Hit &h) {
    // Optional: Change this brute force method into a faster one.
    if (!Box::is_intersect(min, max, r, tmin)) return 0;
    bool result = false;
    for (int triId = 0; triId < (int) t.size(); ++triId) {
        TriangleIndex &triIndex = t[triId];
        Triangle triangle(v[triIndex[0]],
                          v[triIndex[1]], v[triIndex[2]], refl, e, color);
        triangle.normal = n[triId];
        result |= (triangle.intersect(r, tmin, h) != 0);
    }
    return h.t;
}

Mesh::Mesh(const char *filename, Vector3f center, Refl_t refl_, Vector3f e_, Vector3f c_) : Object3D(refl_, e_, c_) {

    // Optional: Use tiny obj loader to replace this simple one.
    std::ifstream f;
    f.open(filename);
    if (!f.is_open()) {
        std::cout << "Cannot open " << filename << "\n";
        return;
    }
    std::string line;
    std::string vTok("v");
    std::string fTok("f");
    std::string texTok("vt");
    char bslash = '/', space = ' ';
    std::string tok;
    int texID;
    while (true) {
        std::getline(f, line);
        if (f.eof()) {
            break;
        }
        if (line.size() < 3) {
            continue;
        }
        if (line.at(0) == '#') {
            continue;
        }
        std::stringstream ss(line);
        ss >> tok;
        if (tok == vTok) {
            Vector3f vec;
            ss >> vec.x() >> vec.y() >> vec.z();
            vec.x() = vec.x() * 100 + center.x();
            vec.y() = vec.y() * 100 + center.y();
            vec.z() = vec.z() * 100 + center.z();
            v.push_back(vec);
        } else if (tok == fTok) {
            if (line.find(bslash) != std::string::npos) {
                std::replace(line.begin(), line.end(), bslash, space);
                std::stringstream facess(line);
                TriangleIndex trig;
                facess >> tok;
                for (int ii = 0; ii < 3; ii++) {
                    facess >> trig[ii] >> texID;
                    trig[ii]--;
                }
                t.push_back(trig);
            } else {
                TriangleIndex trig;
                for (int ii = 0; ii < 3; ii++) {
                    ss >> trig[ii];
                    trig[ii]--;
                }
                t.push_back(trig);
            }
        }
    }
    computeNormal();
    f.close();

    // Bounding Box Intersection Acceleration
    min = {1e6, 1e6, 1e6};
    max = {0, 0, 0};
    for (const auto &i: v) {
        min = {min.x() < i.x() ? min.x() : i.x(), min.y() < i.y() ? min.y() : i.y(), min.z() < i.z() ? min.z() : i.z()};
        max = {max.x() > i.x() ? max.x() : i.x(), max.y() > i.y() ? max.y() : i.y(), max.z() > i.z() ? max.z() : i.z()};
    }
}

void Mesh::computeNormal() {
    n.resize(t.size());
    for (int triId = 0; triId < (int) t.size(); ++triId) {
        TriangleIndex &triIndex = t[triId];
        Vector3f a = v[triIndex[1]] - v[triIndex[0]];
        Vector3f b = v[triIndex[2]] - v[triIndex[0]];
        b = Vector3f::cross(a, b);
        n[triId] = b / ::sqrt(Vector3f::dot(b,b));
    }
}

