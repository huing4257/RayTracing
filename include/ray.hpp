#ifndef RAY_H
#define RAY_H

#include <cassert>
#include <iostream>
#include "vector.h"


// Ray class mostly copied from Peter Shirley and Keith Morley
struct Ray {
    Vec origin, direction;
    Ray(Vec o_, Vec d_) : origin(o_), direction(d_) {}
};


#endif // RAY_H
