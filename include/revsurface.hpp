#ifndef REVSURFACE_HPP
#define REVSURFACE_HPP

#include "object3d.hpp"
#include "curve.hpp"
#include "box.h"
#include <tuple>

class RevSurface : public Object3D {

    Curve *pCurve;

public:
    RevSurface(Curve *pCurve, double x, double z, Refl_t refl_, Vec e_, Vec c_) : pCurve(pCurve), x(x), z(z),
                                                                                  Object3D(refl_, e_, c_) {
        // x always less than 0.
        double min_x = 0, max_y = 1e-6, min_y = 1e6;
        for (const auto &cp: pCurve->getControls()) {
            min_x = std::min(min_x, cp.x);
            max_y = std::max(max_y, cp.y);
            min_y = std::min(min_y, cp.y);
            if (cp.z != 0.0) {
                printf("Profile of revSurface must be flat on xy plane.\n");
                exit(0);
            }
        }
        vmin = Vec(x + min_x, min_y, z + min_x);
        vmax = Vec(x - min_x, max_y, z - min_x);
    }

    ~RevSurface() override {
        delete pCurve;
    }

    double intersect(const Ray &r, double tmin, Hit &h) override {
        double tr;
        if (!Box::is_intersect(vmin, vmax, r, tmin, &tr)) return 0;
        // newton's method
        Vec p = r.origin + r.direction * tr;
        double estimate_t = pCurve->estimate(p.y);
        Vec o = r.origin - Vec(x, 0, z);
        Vec d = r.direction;
        for (int i = 0; i < 10; i++) {
            double f = pCurve->F(o, d, estimate_t);
            double df = pCurve->dF(o, d, estimate_t);
            estimate_t = estimate_t - f / df;
            if (estimate_t < tmin || estimate_t > 1) {
                return 0;
            }
            if (fabs(f) < 1e-6) {
                CurvePoint cp = pCurve->get_pos(estimate_t);
                h.t = (cp.V.y - o.y) / d.y;
                Vec hit_point = o + d * h.t;
                double dis = sqrt(hit_point.x * hit_point.x + hit_point.y * hit_point.y);
                Vec rad = {cp.T.x * hit_point.x / dis, 0, cp.T.x * hit_point.z / dis};
                Vec norm = cp.T % rad % cp.T;
                h.hit_color = color;
                h.refl = refl;
                h.hit_normal = norm;
                return h.t;
            }
        }
        return 0;
    }

    double x, z; // rotation axis
    Vec vmin, vmax;
};

#endif //REVSURFACE_HPP
