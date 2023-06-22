//
// Created by 胡 on 2023/6/21.
//

#ifndef RT_HIT_H
#define RT_HIT_H
#include "vecmath.h"
enum Refl_t {
    DIFF, SPEC, REFR
};  // material types, used in radiance()

struct Hit{
    double t;
    Vector3f hit_pos;
    Vector3f hit_normal;
    Vector3f hit_color;
    Refl_t refl;
    Hit(double t_, Vector3f hitPos_, Vector3f hitNormal_, Vector3f hitColor_, Refl_t refl_):t(t_), hit_pos(hitPos_), hit_normal(hitNormal_), hit_color(hitColor_), refl(refl_){}
    Hit():t(MAXFLOAT), hit_pos({}), hit_normal({}), hit_color({}), refl(DIFF){}
};

Vector3f min(const Vector3f &a, const Vector3f &b){
    return {std::min(a.x(), b.x()), std::min(a.y(), b.y()), std::min(a.z(), b.z())};
}
Vector3f max(const Vector3f &a, const Vector3f &b){
    return {std::max(a.x(), b.x()), std::max(a.y(), b.y()), std::max(a.z(), b.z())};
}
#endif //RT_HIT_H
