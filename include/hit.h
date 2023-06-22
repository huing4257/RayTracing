//
// Created by 胡 on 2023/6/21.
//

#ifndef RT_HIT_H
#define RT_HIT_H
#include <vecmath.h>
enum Refl_t {
    DIFF, SPEC, REFR
};  // material types, used in radiance()

struct Hit{
    double t;
    Vector3f hit_pos;
    Vector3f hit_normal;
    Vector3f hit_color;
    Refl_t refl;
    Hit(double t_, Vector3f hitPos_, Vector3f hitNormal_, Vector3f hitColor_, Refl_t refl_): t(t_), hit_pos(hitPos_), hit_normal(hitNormal_), hit_color(hitColor_), refl(refl_){}
    Hit():t(MAXFLOAT), hit_pos({}), hit_normal({}), hit_color({}), refl(DIFF){}
};
#endif //RT_HIT_H
