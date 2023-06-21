#ifndef PLANE_H
#define PLANE_H

#include "object3d.hpp"
#include "vector.h"
#include <cmath>
#include "texture.h"

// function: ax+by+cz=d
// choose your representation , add more fields and fill in the functions

class Plane : public Object3D {
public:
    Plane() {

    }

    Plane(const Vec &normal, float d, Vec e_, Vec c_, Refl_t refl_,Texture* t = nullptr) : Object3D(refl_,e_,c_) {
        this->normal = normal;
        this->d = d;
        this->texture = t;
    }

    ~Plane() override = default;

    double intersect(const Ray &r, double tmin, Hit &hit) override {// returns distance, 0 if nohit
        float t = (d - Vec::dot(normal, r.origin)) / Vec::dot(normal, r.direction);
        if (t > tmin && t < hit.t) {
            hit.t = t;
            hit.hit_pos = r.origin + r.direction * t;
            hit.hit_normal = normal;
            if(texture != nullptr){
                hit.hit_color = texture->get_color(abs(hit.hit_pos.x), abs(hit.hit_pos.z));
            }
            else { hit.hit_color = color; }
            hit.refl = refl;
            return t;
        }
        return 0;
    }

    Vec normal;
    float d;
    Texture *texture;
protected:
};

#endif //PLANE_H
		

