#ifndef PLANE_H
#define PLANE_H

#include "object3d.hpp"
#include <cmath>
#include "texture.h"
#include <vector>

using std::vector;
// function: ax+by+cz=d
// choose your representation , add more fields and fill in the functions

class Plane : public Object3D {
public:
    Plane() {

    }

    Plane(const Vector3f &normal, float d, Vector3f e_, Vector3f c_, Refl_t refl_, vector<Texture *> t = {}) : Object3D(
            refl_, e_, c_) {
        this->normal = normal;
        this->d = d;
        this->textures = t;
    }

    ~Plane() override = default;

    double intersect(const Ray &r, double tmin, Hit &hit) override {// returns distance, 0 if nohit
        float t = (d - Vector3f::dot(normal, r.origin)) / Vector3f::dot(normal, r.direction);
        if (t > tmin && t < hit.t) {
            hit.t = t;
            hit.hit_pos = r.origin + r.direction * t;
            hit.hit_normal = normal;
            hit.hit_color = color;
            for(auto texture:textures){
                if (texture != nullptr) {
                    if (normal[0] == 0 && normal[2] == 0) {
                        texture->change_hit(abs(hit.hit_pos[0]) / 100, abs(hit.hit_pos[2]) / 170, hit);
                    } else if (normal[0] == 0 && normal[1] == 0) {
                        texture->change_hit(abs(hit.hit_pos[0]) / 100, abs(hit.hit_pos[1]) / 81.6, hit);
                    } else if (normal[1] == 0 && normal[2] == 0) {
                        texture->change_hit(abs(hit.hit_pos[2]) / 170, abs(hit.hit_pos[1]) / 81.6, hit);
                    }
                }
            }
            hit.refl = refl;
            return t;
        }
        return 0;
    }

    Vector3f normal;
    float d;
    vector<Texture *> textures;
protected:
};

#endif //PLANE_H
		

