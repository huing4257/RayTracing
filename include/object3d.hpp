#ifndef OBJECT3D_H
#define OBJECT3D_H

#include "ray.hpp"
#include "hit.h"


// Base class for all 3d entities.
class Object3D {
public:
    Object3D() : refl(DIFF), e({}), color({}) {}

    virtual ~Object3D() = default;

    explicit Object3D(Refl_t refl_, Vector3f e_, Vector3f c_) : refl(refl_), e(e_), color(c_) {}

    // Intersect Ray with this object. If hit, store information in hit structure.
    virtual double intersect(const Ray &r, double tmin, Hit &hit) = 0;

    Refl_t refl;
    Vector3f e; // emission
    Vector3f color; // color
};

#endif

