#ifndef SPHERE_H
#define SPHERE_H

#include "vector.h"
#include "ray.hpp"
#include "object3d.hpp"

struct Sphere : public Object3D {
    double rad; // radius
    Vec center;// position, emission, color
    Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_) : Object3D(refl_, e_,c_), rad(rad_), center(p_) {};

    double intersect(const Ray &r, double tmin, Hit &hit) override {// returns distance, 0 if nohit
        Vec op = center -
                 r.origin;                 // Solve t^2*direction.direction + 2*t*(origin-center).direction + (origin-center).(origin-center)-R^2 = 0
        double t, eps = tmin, b = op.dot(r.direction), det = b * b - op.dot(op) + rad * rad;
        if (det < 0) return 0;
        else
            det = sqrt(det);
        t = b - det;
        if (t > eps && t < hit.t) {
            hit.t = t;
            hit.hit_pos = r.origin + r.direction * t;
            hit.hit_normal = (hit.hit_pos - center).norm();
            hit.hit_color = color;
            hit.refl = refl;
            return t;
        } else {
            t = b + det;
            if (t > eps && t < hit.t) {
                hit.t = t;
                hit.hit_pos = r.origin + r.direction * t;
                hit.hit_normal = (hit.hit_pos - center).norm();
                hit.hit_color = color;
                hit.refl = refl;
                return t;
            }
        }
        return 0;
    }

};

#endif
