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
        double max_x=1e-6, max_y=1e-6;
        // Check flat.
        for (const auto &cp: pCurve->getControls()) {
            max_x = std::max(max_x, cp.x);
            max_y = std::max(max_y, cp.y);
            if (cp.z != 0.0) {
                printf("Profile of revSurface must be flat on xy plane.\n");
                exit(0);
            }
        }
        vmin = Vec(x-max_x, max_y, z-max_x);
        vmax = Vec(x+max_x, max_y, z+max_x);

    }

    ~RevSurface() override {
        delete pCurve;
    }

    double intersect(const Ray &r, double tmin, Hit &h) override {
        if (!Box::is_intersect(vmin, vmax,r,tmin)) return 0;
        double t = solve_t(r, tmin, h);
        if (t <= 0 || t >= 1) return 0;
        CurvePoint pos = pCurve->get_pos(t);
        double tr = (pos.V.y - r.origin.y) / r.direction.y;
        Vec hit_pos = r.origin + r.direction * tr;
        Vec rad = Vec(hit_pos.x - x, 0, hit_pos.z - z).norm();
        Vec new_T = Vec(pos.T.x * rad.x, pos.T.y, pos.T.x * rad.z).norm();
        Vec normal = (new_T % rad) % (new_T).norm();
        if (tr < h.t) {
            h.t = tr;
            h.hit_pos = hit_pos;
            h.hit_normal = normal;
            h.hit_color = color;
            h.refl = refl;
        }
        return t;
    }

    double solve_t(const Ray &r, double tmin, Hit &h) {
        double t = 0.5, ft, dft;

        double ox = r.origin.x - x;
        double oz = r.origin.z - z;
        double oy = r.origin.y;
        double dx = r.direction.x;
        double dz = r.direction.z;
        double dy = r.direction.y;

        for (int i = 10; i--;) {
            if (t < 0) t = 0;
            else if (t > 1) t = 1;
            if (ft < 1e-5) return t;
            CurvePoint pos = pCurve->get_pos(t);
            double ty = (pos.V.y - oy) / dy;
            double xt = pos.V.x;
            double dxt = pos.T.x;
            double dyt = pos.T.y;
            ft = ty * ty * (dx * dx + dz * dz) + 2 * ty * (ox + oz) + ox * ox + oz * oz - xt * xt;
            dft = 2 * ty * (dx * dx + dz * dz) * dyt / dy + 2 * (ox + oz) * dyt / dy - 2 * xt * dxt;
            t -= ft / dft;
        }
        return -1;
    }

    double x, z; // rotation axis
    Vec vmin, vmax;
};

#endif //REVSURFACE_HPP
