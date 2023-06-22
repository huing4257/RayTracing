#ifndef REVSURFACE_HPP
#define REVSURFACE_HPP

#include "object3d.hpp"
#include "curve.hpp"
#include "box.h"
#include <tuple>

double F_(const Vector3f &dir, const Vector3f &origin, double x_t, double y_t) {
    double x1 = dir[1] * dir[1] * x_t * x_t;
    double y1 = y_t - origin[1];
    double x2 = pow((y1 * dir[0] + dir[1] * origin[0]), 2);
    double x3 = pow((y1 * dir[2] + dir[1] * origin[2]), 2);
    return x2 + x3 - x1;
}

double F_grad(const Vector3f &dir, const Vector3f &origin, double x_t, double y_t, double x_grad_t, double y_grad_t) {
    double x1 = 2 * dir[1] * dir[1] * x_t * x_grad_t;
    double y1 = y_t - origin[1];
    double x2 = 2 * dir[0] * y_grad_t * (y1 * dir[0] + dir[1] * origin[0]);
    double x3 = 2 * dir[2] * y_grad_t * (y1 * dir[2] + dir[1] * origin[2]);
    return x2 + x3 - x1;
}

class RevSurface : public Object3D {

    Curve *pCurve;

public:
    RevSurface(Curve *pCurve, double x, double z, Refl_t refl_, Vector3f e_, Vector3f c_) : pCurve(pCurve), x(x), z(z),
                                                                                            Object3D(refl_, e_, c_) {
        // x always less than 0.
        float min_x = 0, max_y = 1e-6, min_y = 1e6;
        for (const auto &cp: pCurve->getControls()) {
            min_x = std::min(min_x, cp[0]);
            max_y = std::max(max_y, cp[1]);
            min_y = std::min(min_y, cp[1]);
            if (cp[2] != 0.0) {
                printf("Profile of revSurface must be flat on xy plane.\n");
                exit(0);
            }
        }
        vmin = Vector3f(x + min_x, min_y, z + min_x);
        vmax = Vector3f(x - min_x, max_y, z - min_x);
    }

    ~RevSurface() override {
        delete pCurve;
    }

    double intersect(const Ray &r, double tmin, Hit &h) override {
        double tr;
        if (!Box::is_intersect(vmin, vmax, r, tmin, &tr)) return 0;
        // newton's method
        Vector3f p = r.origin + r.direction * tr;
        double estimate_t = pCurve->estimate(p[1]);
        Vector3f o = r.origin - Vector3f(x, 0, z);
        Vector3f d = r.direction;
        for (int i = 0; i < 10; i++) {
            Vector3f now_point;
            Vector3f now_grad;
            CurvePoint cp = pCurve->get_pos(estimate_t);
            now_point = cp.V;
            now_grad = cp.T;
            double ft = F_(d, o, now_point[0], now_point[1]);
            double ft_grad = F_grad(r.direction, r.origin, now_point[0], now_point[1], now_grad[0], now_grad[1]);
            // printf("dep: %d", dep);
            // printf("ft: %f\n", ft);
            // printf("ft_grad: %f\n", ft_grad);
            // printf("estimate: %f\n", estimate_t);

            if (abs(ft) < 1e-5) {
                if (abs(r.direction[1]) > 1e-3) {
                    tr = (now_point[1] - r.origin[1]) / (r.direction[1]);
                } else {

                    double a = r.direction[2] * r.direction[2] + r.direction[0] + r.direction[0];
                    double b = 2 * (r.direction[2] * r.origin[2] + r.direction[0] * r.origin[0]);
                    double c = pow(r.origin[0], 2) + pow(r.origin[2], 2) - pow(now_point[0], 2);
                    tr = (sqrt(pow(b, 2) - 4 * a * c) - b) / (2 * a);
                }
                Vector3f next_origin = r.direction * tr + r.origin;

                if (tr > tmin && pCurve->is_on_curve(next_origin)) {

                    Vector3f plane_normal(now_grad[1], -now_grad[0], 0);
                    plane_normal.normalized();

                    if (abs(now_point[0]) < 1e-4) {
                        h.t = tr;
                        h.hit_normal = now_point[1] > 0 ? Vector3f(0, 1, 0) : Vector3f(0, -1, 0);
                        h.refl = refl;
                        h.hit_color = color;
                        h.hit_pos = r.origin + r.direction * tr;
                        return tr;
                    }
                    double costheta = (r.direction[0] * tr + r.origin[0]) / now_point[0];
                    double sintheta = (r.direction[2] * tr + r.origin[2]) / now_point[0];
                    Vector3f new_normal = Vector3f(costheta * plane_normal[0], plane_normal[1], sintheta * plane_normal[0]);
                    h.t = tr;
                    h.refl = refl;
                    h.hit_color = color;
                    h.hit_pos = r.origin + r.direction * tr;
                    return tr;
                }
            }
            double step = ft / ft_grad;
            if (step > 0.05)
                step = 0.05;
            else if (step < -0.05)
                step = -0.05;
            estimate_t -= step;
        }
        return 0;
    }

    double x, z; // rotation axis
    Vector3f vmin, vmax;
};

#endif //REVSURFACE_HPP
