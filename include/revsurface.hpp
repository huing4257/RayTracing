#ifndef REVSURFACE_HPP
#define REVSURFACE_HPP

#include "object3d.hpp"
#include "curve.hpp"
#include "box.h"
#include <tuple>

double F_(const Vec &dir, const Vec &origin, double x_t, double y_t) {
    double x1 = dir.y * dir.y * x_t * x_t;
    double y1 = y_t - origin.y;
    double x2 = pow((y1 * dir.x + dir.y * origin.x), 2);
    double x3 = pow((y1 * dir.z + dir.y * origin.z), 2);
    return x2 + x3 - x1;
}

double F_grad(const Vec &dir, const Vec &origin, double x_t, double y_t, double x_grad_t, double y_grad_t) {
    double x1 = 2 * dir.y * dir.y * x_t * x_grad_t;
    double y1 = y_t - origin.y;
    double x2 = 2 * dir.x * y_grad_t * (y1 * dir.x + dir.y * origin.x);
    double x3 = 2 * dir.z * y_grad_t * (y1 * dir.z + dir.y * origin.z);
    return x2 + x3 - x1;
}

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
            Vec now_point;
            Vec now_grad;
            CurvePoint cp = pCurve->get_pos(estimate_t);
            now_point = cp.V;
            now_grad = cp.T;
            double ft = F_(d, o, now_point.x, now_point.y);
            double ft_grad = F_grad(r.direction, r.origin, now_point.x, now_point.y, now_grad.x, now_grad.y);
            // printf("dep: %d", dep);
            // printf("ft: %f\n", ft);
            // printf("ft_grad: %f\n", ft_grad);
            // printf("estimate: %f\n", estimate_t);

            if (abs(ft) < 1e-5) {
                double tr;
                if (abs(r.direction.y) > 1e-3) {
                    tr = (now_point.y - r.origin.y) / (r.direction.y);
                    // if(tr < 0) {

                    // }
                } else {

                    double a = r.direction.z * r.direction.z + r.direction.x + r.direction.x;
                    double b = 2 * (r.direction.z * r.origin.z + r.direction.x * r.origin.x);
                    double c = pow(r.origin.x, 2) + pow(r.origin.z, 2) - pow(now_point.x, 2);
                    tr = (sqrt(pow(b, 2) - 4 * a * c) - b) / (2 * a);
                    // std::cout << tr << " " << a  << " "<< b << " " << c << std::endl;
                    // now_point.print();
                    // r.origin.print();
                    // r.direction.print();
                }
                Vec next_origin = r.direction * tr + r.origin;

                if (tr > tmin && pCurve->is_on_curve(next_origin)) {
                    //  std::cout << AABB_x << " " << AABB_y << std::endl;
                    //     now_point.print();
                    //     r.origin.print();
                    //     r.direction.print();


                    Vec plane_normal(now_grad.y, -now_grad.x, 0);
                    plane_normal.norm();

                    if (abs(now_point.x) < 1e-4) {
                        h.t = tr;
                        h.hit_normal = now_point.y > 0 ? Vec(0, 1, 0) : Vec(0, -1, 0);
                        h.refl = refl;
                        h.hit_color = color;
                        h.hit_pos = next_origin;
                        return tr;
                    }
                    double costheta = (r.direction.x * tr + r.origin.x) / now_point.x;
                    double sintheta = (r.direction.z * tr + r.origin.z) / now_point.x;
                    Vec new_normal = Vec(costheta * plane_normal.x, plane_normal.y, sintheta * plane_normal.x);
                    // if(Vec::dot(new_normal, next_origin) < 0)
                    //     new_normal.negate();
                    // std::cout << "nextori";
                    // next_origin.print();
                    // std::cout << "normal";
                    // new_normal.print();
                    h.t = tr;
                    h.refl = refl;
                    h.hit_color = color;
                    return true;
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
    Vec vmin, vmax;
};

#endif //REVSURFACE_HPP
