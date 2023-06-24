#ifndef REVSURFACE_HPP
#define REVSURFACE_HPP

#include "object3d.hpp"
#include "curve.hpp"
#include "box.h"
#include <tuple>
#include <float.h>

const int resolution = 10;
const int NEWTON_STEPS = 20;
const float NEWTON_EPS = 1e-4;

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
        vmin = Vector3f(  min_x, min_y, min_x);
        vmax = Vector3f( - min_x, max_y,- min_x);
    }

    ~RevSurface() override {
        delete pCurve;
    }

    double intersect(const Ray &r, double tmin, Hit &hit) override {
        // return meshIntersect(r, hit, tmin);
        Ray new_ray = r;
        new_ray.origin -= Vector3f(x,0,z);
        if (!newtonIntersect(new_ray, hit)) return 0;
        hit.hit_pos = r.origin + hit.t*r.direction;
        return hit.t;
    }

    bool newtonIntersect(const Ray &r, Hit &hit) {
        float t , theta, mu;
        if (!Box::is_intersect(vmin, vmax, r, 1e-5,&t)||t > hit.t) return false;
        getUV(r, t, theta, mu);
        Vector3f normal, point;
        if (!newton(r, t, theta, mu, normal, point)) {
            return false;
        }
        if (!isnormal(mu) || !isnormal(theta) || !isnormal(t)) return false;
        if (t < 0 || mu < pCurve->range[0] || mu > pCurve->range[1] ||
            t > hit.t)
            return false;
        hit = {t, point, normal.normalized(), color, refl};
        return true;
    }

    bool newton(const Ray &r, float &t, float &theta, float &mu,
                Vector3f &normal, Vector3f &point) {
        Vector3f dmu, dtheta;
        for (int i = 0; i < 10; ++i) {
            if (theta < 0.0) theta += 2 * M_PI;
            if (theta >= 2 * M_PI) theta = fmod(theta, 2 * M_PI);
            if (mu >= 1) mu = 1.0 - FLT_EPSILON;
            if (mu <= 0) mu = FLT_EPSILON;
            point = getPoint(theta, mu, dtheta, dmu);
            Vector3f f = r.origin + r.direction * t - point;
            float dist2 = f.squaredLength();
            normal = Vector3f::cross(dmu, dtheta);
            if (dist2 < NEWTON_EPS) return true;
            float D = Vector3f::dot(r.direction, normal);
            t -= Vector3f::dot(dmu, Vector3f::cross(dtheta, f)) / D;
            mu -= Vector3f::dot(r.direction, Vector3f::cross(dtheta, f)) / D;
            theta += Vector3f::dot(r.direction, Vector3f::cross(dmu, f)) / D;
        }
        return false;
    }

    Vector3f getPoint(const float &theta, const float &mu, Vector3f &dtheta,
                      Vector3f &dmu) {
        Vector3f pt;
        Quat4f rot;
        rot.setAxisAngle(theta, Vector3f::UP);
        Matrix3f rotMat = Matrix3f::rotation(rot);
        CurvePoint cp = pCurve->getPoint(mu);
        pt = rotMat * cp.V;
        dmu = rotMat * cp.T;
        dtheta = Vector3f(-cp.V.x() * sin(theta), 0, -cp.V.x() * cos(theta));
        return pt;
    }

    void getUV(const Ray &r, const float &t, float &theta, float &mu) {
        Vector3f pt(r.origin + r.direction * t);
        theta = atan2(-pt.z(), pt.x()) + M_PI;
        mu = (pCurve->ymax - pt.y()) / (pCurve->ymax - pCurve->ymin);
    }

    double x, z; // rotation axis
    Vector3f vmin, vmax;
};

#endif //REVSURFACE_HPP
