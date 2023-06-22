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
};

class Curve : public Object3D {
protected:
    std::vector<Vec> controls;
public:
    explicit Curve(std::vector<Vec> points) : controls(std::move(points)) {}

    double intersect(const Ray &r, double tmin, Hit &h) override {
        return false;
    }

    std::vector<Vec> &getControls() {
        return controls;
    }

    virtual CurvePoint get_pos(double t) = 0;

};

class BezierCurve : public Curve {
public:
    explicit BezierCurve(const std::vector<Vec> &points) : Curve(points) {
        if (points.size() < 4 || points.size() % 3 != 1) {
            printf("Number of control points of BezierCurve must be 3n+1!\n");
            exit(0);
        }
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

class BsplineCurve : public Curve {
public:
    BsplineCurve(const std::vector<Vec> &points) : Curve(points) {
        if (points.size() < 4) {
            printf("Number of control points of BspineCurve must be more than 4!\n");
            exit(0);
        }
        n = points.size() - 1;
        knots.resize(n + k + 2);
        for (int i = 0; i <= n + k + 1; ++i) {
            knots[i] = (double) i / (double) (n + k + 1);
        }
    }

    CurvePoint get_pos(double t) override {
        auto dp = new double *[n + k + 2];
        for (int l = 0; l < n + k + 1; ++l) {
            dp[l] = new double[k];
            for (int p = 0; p < k; ++p) {
                dp[l][p] = 0;
            }
        }
        for (int l = 0; l < n + k + 1; ++l) {
            dp[l][0] = get_B_i_0(l, t);
        }
        for (int p = 1; p <= k; p++) {
            for (int l = 0; l < n + k + 2; ++l) {
                if (l + p > n + k) {
                    continue;
                }
                dp[l][p] = (t - knots[l]) / (knots[l + p] - knots[l]) * dp[l][p - 1] +
                           (knots[l + p + 1] - t) / (knots[l + p + 1] - knots[l + 1]) * dp[l + 1][p - 1];
            }
        }

        CurvePoint cp;
        cp.V = Vec(0, 0, 0);
        for (int l = 0; l <= n; ++l) {
            cp.V += controls[l] * dp[l][k];
        }
        cp.T = Vec(0, 0, 0);
        for (int l = 0; l <= n; ++l) {
            cp.T += controls[l] * (double) k * ((dp[l][k - 1] / (knots[l + k] - knots[l])) -
                                                (dp[l + 1][k - 1] / (knots[l + k + 1] - knots[l + 1])));
        }
        cp.T.norm();
        delete[] dp;
        return cp;
    }

protected:
    const int k = 3;
    int n;
    std::vector<double> knots;

    double get_B_i_0(int i, double t) {
        if (t >= knots[i] && t < knots[i + 1]) {
            return 1;
        } else {
            return 0;
        }
    }
};

#endif // CURVE_HPP
