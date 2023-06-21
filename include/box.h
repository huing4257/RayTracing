//
// Created by 胡 on 2023/6/21.
//

#ifndef RT_BOX_H
#define RT_BOX_H

#include "object3d.hpp"

class Box{
public:
    static bool is_intersect(Vec min, Vec max,const Ray &r, double tmin, Hit &h) {
        double t11 = (min.x - r.origin.x) / r.direction.x;
        double t12 = (max.x - r.origin.x) / r.direction.x;
        double t13 = (min.y - r.origin.y) / r.direction.y;
        double t21 = (max.y - r.origin.y) / r.direction.y;
        double t22 = (min.z - r.origin.z) / r.direction.z;
        double t23 = (max.z - r.origin.z) / r.direction.z;
        double t1 = std::min(t11, t12);
        t1 = std::min(t1, t13);
        double t2 = std::max(t21, t22);
        t2 = std::max(t2, t23);
        if(t1 >= t2 || t2 < tmin) return false;
        return true;
    }
};

#endif //RT_BOX_H
