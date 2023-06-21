//
// Created by 胡 on 2023/6/21.
//

#ifndef RT_VECTOR_H
#define RT_VECTOR_H

struct Vec {        // Usage: time ./smallpt4k && xv image.ppm
    double x, y, z;                  // position, also color (r,g,b)
    Vec(double x_ = 0, double y_ = 0, double z_ = 0) {
        x = x_;
        y = y_;
        z = z_;
    }

    Vec operator+(const Vec &b) const { return {x + b.x, y + b.y, z + b.z}; }

    Vec operator-(const Vec &b) const { return {x - b.x, y - b.y, z - b.z}; }

    Vec operator*(double b) const { return {x * b, y * b, z * b}; }

    Vec operator/(double b) const { return {x / b, y / b, z / b}; }

    Vec mult(const Vec &b) const { return {x * b.x, y * b.y, z * b.z}; }

    Vec &norm() { return *this = *this * (1 / sqrt(x * x + y * y + z * z)); }

    double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; }

    // cross:
    Vec operator%(Vec &b) const { return {y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x}; }

    static double dot(const Vec &a, const Vec &b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

    static Vec cross(const Vec &a, const Vec &b) {
        return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
    }

};

#endif //RT_VECTOR_H
