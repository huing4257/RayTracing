#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>
#include <iostream>
using namespace std;

// TODO: implement this class and add more fields as necessary,
class Triangle: public Object3D {

public:
	Triangle() = delete;

    // a b c are three vertex positions of the triangle
    Triangle(const Vector3f &a, const Vector3f &b, const Vector3f &c, Material m) : Object3D(m) {
        vertices[0] = a;
        vertices[1] = b;
        vertices[2] = c;
        normal = Vector3f::cross(b - a, c - a).normalized();
        material = m;
    }

    bool intersect(const Ray &ray, Hit &hit, float tmin) override {
        Vector3f e1 = vertices[0] - vertices[1];
        Vector3f e2 = vertices[0] - vertices[2];
        Vector3f s = vertices[0] - ray.getOrigin();
        float det = Vector3f::dot((Vector3f::cross(ray.getDirection(), e1)), e2);
        float t = Vector3f::dot((Vector3f::cross(s, e1)), e2) / det;
        float beta = Vector3f::dot((Vector3f::cross(ray.getDirection(), e1)), s) / det;
        float gamma = Vector3f::dot((Vector3f::cross(s, e2)), ray.getDirection()) / det;
        float alpha = 1 - beta - gamma;
        if (t > tmin && t < hit.getT() && alpha >= 0 && beta >= 0 && gamma >= 0) {
            hit.set(t, material, normal);
            return true;
        }
        return false;
	}
	Vector3f normal;
	Vector3f vertices[3];
protected:

};

#endif //TRIANGLE_H
