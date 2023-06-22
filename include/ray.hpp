#ifndef RAY_H
#define RAY_H

#include <cassert>
#include <iostream>
#include <vecmath.h>


// Ray class mostly copied from Peter Shirley and Keith Morley
struct Ray {
    Vector3f origin, direction;
    Ray(Vector3f o_, Vector3f d_) : origin(o_), direction(d_) {}
};


#endif // RAY_H
