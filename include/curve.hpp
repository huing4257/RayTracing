#ifndef CURVE_HPP
#define CURVE_HPP

#include "object3d.hpp"
#include "vector.h"
#include <vector>
#include <utility>

#include <algorithm>

// TODO (PA2): Implement Bernstein class to compute spline basis function.
//       You may refer to the python-script for implementation.

// The CurvePoint object stores information about a point on a curve
// after it has been tesselated: the vertex (V) and the tangent (T)
// It is the responsiblility of functions that create these objects to fill in all the data.
struct CurvePoint {
    Vec V; // Vertex
    Vec T; // Tangent  (not normalized)
    double t;
};

class Curve : public Object3D {
protected:
    std::vector<Vec> controls;
    std::vector<CurvePoint> data;
public:
    explicit Curve(std::vector<Vec> points) : controls(std::move(points)) {}

    double intersect(const Ray &r, double tmin, Hit &h) override {
        return false;
    }

    virtual double estimate(double y) = 0;

    std::vector<Vec> &getControls() {
        return controls;
    }

    virtual CurvePoint get_pos(double t) = 0;

    virtual double F(Vec &o, Vec &d, double t) = 0;

    virtual double dF(Vec &o, Vec &d, double t) = 0;
};

class BezierCurve : public Curve {
public:
    explicit BezierCurve(const std::vector<Vec> &points) : Curve(points) {
        if (points.size() < 4 || points.size() % 3 != 1) {
            printf("Number of control points of BezierCurve must be 3n+1!\n");
            exit(0);
        }
        for (int i = 0; i < 100; ++i) {
            float t = (float) i / (float) 100;
            size_t n = controls.size();
            auto *mid_points = new CurvePoint[n];
            for (int j = 0; j < n; ++j) {
                mid_points[j].V = controls[j];
            }
            for (size_t k = n - 1; k > 0; k--) {
                if (k == 1) {
                    mid_points[0].T = mid_points[1].V - mid_points[0].V;
                }
                for (int j = 0; j < k; ++j) {
                    mid_points[j].V = mid_points[j].V * t + mid_points[j + 1].V * (1 - t);
                }
            }
            mid_points[0].t = t;
            data.push_back(mid_points[0]);
            delete[] mid_points;
        }
    }

    double F(Vec &o, Vec &d, double t) override {
        CurvePoint cp = get_pos(t);
        double t_ = (cp.V.y - o.y) / d.y;
        return t_ * t_ * (d.x * d.x + d.z * d.z) + o.x * o.x + o.z * o.z + 2 * (o.x + o.z) * t_ - cp.V.x * cp.V.x;
    }

    double dF(Vec &o, Vec &d, double t) override {
        CurvePoint cp = get_pos(t);
        double t_ = (cp.V.y - o.y) / d.y;
        return 2 * t_ * (d.x * d.x + d.z * d.z)/d.y*cp.T.y + 2 * (o.x + o.z)/d.y*cp.T.y - 2 * cp.V.x * cp.T.x;
    }



    double estimate(double y) override {
        return min_element(data.begin(), data.end(), [y](const CurvePoint &a, const CurvePoint &b) {
            return abs(a.t - y) < abs(b.t - y);
        })->t;
    }

    CurvePoint get_pos(double t) override {
        size_t n = controls.size();
        auto *mid_points = new CurvePoint[n];
        for (int j = 0; j < n; ++j) {
            mid_points[j].V = controls[j];
        }
        for (size_t k = n - 1; k > 0; k--) {
            if (k == 1) {
                mid_points[0].T = mid_points[1].V - mid_points[0].V;
            }
            for (int j = 0; j < k; ++j) {
                mid_points[j].V = mid_points[j].V * t + mid_points[j + 1].V * (1 - t);
            }
        }
        CurvePoint res = mid_points[0];
        delete[] mid_points;
        return res;
    }

protected:

};

#endif // CURVE_HPP
