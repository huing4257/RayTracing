#ifndef CURVE_HPP
#define CURVE_HPP

#include "object3d.hpp"
#include <vecmath.h>
#include <vector>
#include <utility>

#include <algorithm>

// TODO (PA2): Implement Bernstein class to compute spline basis function.
//       You may refer to the python-script for implementation.

// The CurvePoint object stores information about a point on a curve
// after it has been tesselated: the vertex (V) and the tangent (T)
// It is the responsiblility of functions that create these objects to fill in all the data.
struct CurvePoint {
    Vector3f V; // Vertex
    Vector3f T; // Tangent  (unit)
};

class Curve : public Object3D {
protected:
    std::vector<Vector3f> controls;
public:
    explicit Curve(std::vector<Vector3f> points) : controls(std::move(points)) {}

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        return false;
    }

    std::vector<Vector3f> &getControls() {
        return controls;
    }

    virtual void discretize(int resolution, std::vector<CurvePoint> &data) = 0;

    void drawGL() override {
        Object3D::drawGL();
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glColor3f(1, 1, 0);
        glBegin(GL_LINE_STRIP);
        for (auto &control: controls) { glVertex3fv(control); }
        glEnd();
        glPointSize(4);
        glBegin(GL_POINTS);
        for (auto &control: controls) { glVertex3fv(control); }
        glEnd();
        std::vector<CurvePoint> sampledPoints;
        discretize(30, sampledPoints);
        glColor3f(1, 1, 1);
        glBegin(GL_LINE_STRIP);
        for (auto &cp: sampledPoints) { glVertex3fv(cp.V); }
        glEnd();
        glPopAttrib();
    }
};

class BezierCurve : public Curve {
public:
    explicit BezierCurve(const std::vector<Vector3f> &points) : Curve(points) {
        if (points.size() < 4 || points.size() % 3 != 1) {
            printf("Number of control points of BezierCurve must be 3n+1!\n");
            exit(0);
        }
    }

    void discretize(int resolution, std::vector<CurvePoint> &data) override {
        data.clear();
        // TODO (PA2): fill in data vector
        for (int i = 0; i < resolution; ++i) {
            float t = (float) i / (float) resolution;
            size_t n = controls.size();
            auto *mid_points = new CurvePoint[n];
            for (int j = 0; j < n; ++j) {
                mid_points[j].V = controls[j];
            }
            for (size_t k = n - 1; k > 0; k--) {
                if (k == 1) {
                    mid_points[0].T = mid_points[1].V - mid_points[0].V;
                    mid_points[0].T.normalize();
                }
                for (int j = 0; j < n; ++j) {
                    mid_points[j].V = t * mid_points[j].V + (1 - t) * mid_points[j + 1].V;
                }
            }
            data.push_back(mid_points[0]);
        }
    }

protected:

};

class BsplineCurve : public Curve {
public:
    BsplineCurve(const std::vector<Vector3f> &points) : Curve(points) {
        if (points.size() < 4) {
            printf("Number of control points of BspineCurve must be more than 4!\n");
            exit(0);
        }
        n = points.size() - 1;
        knots.resize(n + k + 2);
        for (int i = 0; i <= n + k + 1; ++i) {
            knots[i] = (float) i / (float) (n + k + 1);
        }
    }

    void discretize(int resolution, std::vector<CurvePoint> &data) override {
        data.clear();
        // TODO (PA2): fill in data vector
        for (int i = k; i <= n; i++) {
            for (int j = 0; j < resolution; j++) {
                float t = knots[i] + (knots[i + 1] - knots[i]) * (float) j / (float) resolution;

                auto dp = new float *[n + k + 2];
                for (int l = 0; l < n + k + 1; ++l) {
                    dp[l] = new float[k];
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
                cp.V = Vector3f(0, 0, 0);
                for (int l = 0; l <= n; ++l) {
                    cp.V += dp[l][k] * controls[l];
                }
                cp.T = Vector3f(0, 0, 0);
                for (int l = 0; l <= n; ++l) {
                    cp.T += (float) k * ((dp[l][k - 1] / (knots[l + k] - knots[l])) -
                                         (dp[l + 1][k - 1] / (knots[l + k + 1] - knots[l + 1]))) *controls[l];
                }
//                cp.T *= (float) (k - 1);
                data.push_back(cp);
            }
        }
    }

protected:
    const int k = 3;
    int n;
    std::vector<float> knots;

    float get_B_i_0(int i, float t) {
        if (t >= knots[i] && t < knots[i + 1]) {
            return 1;
        } else {
            return 0;
        }
    }
};

#endif // CURVE_HPP
