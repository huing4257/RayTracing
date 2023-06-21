#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "object3d.hpp"
#include "vector.h"
#include <cmath>
#include <iostream>

using namespace std;

// TODO: implement this class and add more fields as necessary,
class Triangle : public Object3D {

public:
    Triangle() = delete;

    // a b c are three vertex positions of the triangle
    Triangle(const Vec &a, const Vec &b, const Vec &c, Refl_t refl, Vec e_, Vec color) : Object3D(refl, e_, color) {
        vertices[0] = a;
        vertices[1] = b;
        vertices[2] = c;
        normal = Vec::cross(b - a, c - a).norm();
    }

    double intersect(const Ray &ray, double tmin, Hit &hit) override {
        Vec e1 = vertices[0] - vertices[1];
        Vec e2 = vertices[0] - vertices[2];
        Vec s = vertices[0] - ray.origin;
        float det = Vec::dot((Vec::cross(ray.direction, e1)), e2);
        float t = Vec::dot((Vec::cross(s, e1)), e2) / det;
        float beta = Vec::dot((Vec::cross(ray.direction, e1)), s) / det;
        float gamma = Vec::dot((Vec::cross(s, e2)), ray.direction / det);
        float alpha = 1 - beta - gamma;
        if (t > tmin && t < hit.t && alpha >= 0 && beta >= 0 && gamma >= 0) {
            hit.t = t;
            hit.hit_pos = ray.origin + ray.direction * t;
            hit.hit_normal = normal;
            hit.hit_color = color;
            hit.refl = refl;
            return t;
        }
        return false;
    }

    Vec normal;
    Vec vertices[3];
protected:

};

#endif //TRIANGLE_H
