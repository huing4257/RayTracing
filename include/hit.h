//
// Created by 胡 on 2023/6/21.
//

#ifndef RT_HIT_H
#define RT_HIT_H
#include "vector.h"
enum Refl_t {
    DIFF, SPEC, REFR
};  // material types, used in radiance()

struct Hit{
    double t;
    Vec hit_pos;
    Vec hit_normal;
    Vec hit_color;
    Refl_t refl;
    Hit(double t_, Vec hitPos_, Vec hitNormal_, Vec hitColor_, Refl_t refl_):t(t_), hit_pos(hitPos_), hit_normal(hitNormal_), hit_color(hitColor_), refl(refl_){}
    Hit():t(MAXFLOAT), hit_pos({}), hit_normal({}), hit_color({}), refl(DIFF){}
};
#endif //RT_HIT_H
