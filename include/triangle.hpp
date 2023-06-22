#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "object3d.hpp"
#include "vecmath.h"
#include <cmath>
#include <iostream>

using namespace std;

// TODO: implement this class and add more fields as necessary,
class Triangle : public Object3D {

public:
    Triangle() = delete;

    // a b c are three vertex positions of the triangle
    Triangle(const Vector3f &a, const Vector3f &b, const Vector3f &c, Refl_t refl, Vector3f e_, Vector3f color)
            : Object3D(refl, e_, color) {
        vertices[0] = a;
        vertices[1] = b;
        vertices[2] = c;
        normal = Vector3f::cross(b - a, c - a).normalized();
    }

    double intersect(const Ray &ray, double tmin, Hit &hit) override {
        Vector3f e1 = vertices[0] - vertices[1];
        Vector3f e2 = vertices[0] - vertices[2];
        Vector3f s = vertices[0] - ray.origin;
        float det = Vector3f::dot((Vector3f::cross(ray.direction, e1)), e2);
        float t = Vector3f::dot((Vector3f::cross(s, e1)), e2) / det;
        float beta = Vector3f::dot((Vector3f::cross(ray.direction, e1)), s) / det;
        float gamma = Vector3f::dot((Vector3f::cross(s, e2)), ray.direction / det);
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

    Vector3f normal;
    Vector3f vertices[3];
protected:

};

#endif //TRIANGLE_H
