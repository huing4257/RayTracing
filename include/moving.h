//
// Created by 胡 on 2023/6/23.
//

#ifndef RT_MOVING_H
#define RT_MOVING_H
#include "object3d.hpp"

class Moving :public Object3D {
public:
    Object3D* obj;
    Vector3f v;
    Moving(Object3D* obj_,Vector3f v_) : Object3D(), obj(obj_),v(v_) {}

    double intersect(const Ray &r, double tmin, Hit &hit) override{
        float time = rand()%1000/1000.0;
        Ray r1(r.origin - v * time, r.direction);
        double t= obj->intersect(r1, tmin, hit);
        hit.hit_pos = hit.hit_pos + v * time;
        return t;
    }
};
#endif //RT_MOVING_H
