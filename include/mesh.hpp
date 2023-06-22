#ifndef MESH_H
#define MESH_H

#include <vector>
#include "object3d.hpp"
#include "triangle.hpp"
#include "vector.h"


class Mesh : public Object3D {

public:
    Mesh(const char *filename,Vec center, Refl_t refl_, Vec e_, Vec c_);

    struct TriangleIndex {
        TriangleIndex() {
            x[0] = 0;
            x[1] = 0;
            x[2] = 0;
        }

        int &operator[](const int i) { return x[i]; }

        // By Computer Graphics convention, counterclockwise winding is front face
        int x[3]{};
    };

    std::vector<Vec> v;
    std::vector<TriangleIndex> t;
    std::vector<Vec> n;
    Vec min, max;
    double intersect(const Ray &r, double tmin, Hit &h) override;

private:

    // Normal can be used for light estimation
    void computeNormal();
};

#endif
