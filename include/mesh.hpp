#ifndef MESH_H
#define MESH_H

#include <vector>
#include "object3d.hpp"
#include "triangle.hpp"


class Mesh : public Object3D {

public:
    Mesh(const char *filename, Vector3f center, Refl_t refl_, Vector3f e_, Vector3f c_);

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

    std::vector<Vector3f> v;
    std::vector<TriangleIndex> t;
    std::vector<Vector3f> n;
    Vector3f min, max;
    double intersect(const Ray &r, double tmin, Hit &h) override;

private:

    // Normal can be used for light estimation
    void computeNormal();
};

#endif
