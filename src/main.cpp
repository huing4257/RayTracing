#include <cmath>  // smallpt, a Path Tracer by Kevin Beason, 2008
#include <cstdio> //        Remove "-fopenmp" for g++ version < 4.2
#include <cstdlib>// Make : g++ -O3 -fopenmp smallpt.cpp -origin smallpt
#include "sphere.hpp"
#include "ray.hpp"
#include "object3d.hpp"
#include "hit.h"
#include "plane.hpp"
#include "triangle.hpp"
#include "mesh.hpp"
#include "revsurface.hpp"
#include <iostream>

#define SCALAR 10

BezierCurve pCurve(std::vector<Vector3f>{Vector3f(-0.9, 1.6, 0) * SCALAR,
                                         Vector3f(-1, 1.5, 0) * SCALAR,
                                         Vector3f(-1, 0.8, 0) * SCALAR,
                                         Vector3f(-0.6, 0.4, 0) * SCALAR,
                                         Vector3f(-0.4, 0.4, 0) * SCALAR,
                                         Vector3f(-0.2, 0.15, 0) * SCALAR,
                                         Vector3f(-0.75, 0, 0) * SCALAR,});


Mapping texture = Mapping("../texture/wood.jpg");
NormalMapping normalTexture = NormalMapping("../texture/Wall_n.png");

Object3D *objs[] = {
        new Plane(Vector3f(1, 0, 0), 0, Vector3f(0, 0, 0), Vector3f(.75, .75, .75), DIFF, &normalTexture),  //Left
        new Plane(Vector3f(-1, 0, 0), -100, Vector3f(0, 0, 0), Vector3f(.75, .75, .75), DIFF, &normalTexture),//Rght
        new Plane(Vector3f(0, 0, 1), 0, Vector3f(0,0,0), Vector3f(.75, .75, .75), DIFF),        //Back
        new Plane(Vector3f(0, 0, -1), -170, Vector3f(0,0,0), Vector3f(0,0,0), DIFF, &normalTexture),              //Frnt
        new Plane(Vector3f(0, 1, 0), 0, Vector3f(0,0,0), Vector3f(.75, .75, .75), DIFF, &texture),        //Botm
        new Plane(Vector3f(0, -1, 0), -81.6, Vector3f(0,0,0), Vector3f(.75, .75, .75), DIFF),//Top
        new Sphere(600, Vector3f(50, 681.6 - .27, 81.6), Vector3f(12, 12, 12), Vector3f(), DIFF),   //Lite
        // new Sphere(10, Vector3f(20, 10, 50), Vector3f(0,0,0), Vector3f(.99, .99, .99), REFR),
        // new Mesh("../mesh/bunny_200.obj", Vector3f(70, 5, 80), DIFF, Vector3f(0, 0, 0), Vector3f(.75, .75, .75)),
        new RevSurface(&pCurve, 30, 75, DIFF, Vector3f(0, 0, 0), Vector3f(.75, .75, .75) )
};

inline double clamp(double x) {
    return x < 0 ? 0 : x > 1 ? 1 : x;
}

inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }

inline bool intersect(const Ray &r, double &t, int &id, Hit &hit) {
    double num = (int) (sizeof(objs) / sizeof(Object3D *));
    double d, inf = t = 1e20;
    for (int i = int(num); i--;)
        if ((d = objs[i]->intersect(r, 1e-4, hit)) != 0 && d < t) {
            t = d;
            id = i;
        }
    return t < inf;
}

Vector3f radiance(const Ray &r, int depth, unsigned short *Xi) {
    double t; // distance to intersection
    int id = 0;                            // id of intersected object
    Hit hit;
    if (!intersect(r, t, id, hit)) return {};// if miss, return black
    const Object3D *obj = objs[id];       // the hit object
    Vector3f n = hit.hit_normal, nl = Vector3f::dot(r.direction, n) < 0 ? n : n * -1;
    Vector3f x = hit.hit_pos, f = hit.hit_color;
    double p = f.x() > f.y() && f.x() > f.z() ? f.x() : f.y() > f.z() ? f.y()
                                                                      : f.z();// max refl
    if (++depth > 5) {
        if (erand48(Xi) < p) f = f * (1 / p);
        else
            return obj->e;  //R.R.
    }

    if (obj->refl == DIFF) {
        // Ideal DIFFUSE reflection
        double r1 = 2 * M_PI * erand48(Xi), r2 = erand48(Xi), r2s = sqrt(r2);
        Vector3f w = nl, u = Vector3f::cross((fabs(w.x()) > .1 ? Vector3f(0, 1, 0) : Vector3f(1, 0, 0)),
                                             w).normalized(),
                v = Vector3f::cross(w, u);
        Vector3f d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).normalized();
        return obj->e + f * radiance(Ray(x, d), depth, Xi);

    } else if (obj->refl == SPEC) {// Ideal SPECULAR reflection
        return obj->e + f * (radiance(Ray(x, r.direction - n * 2 * Vector3f::dot(n, r.direction)), depth, Xi));
    }

    Ray reflRay(x, r.direction - n * 2 * Vector3f::dot(r.direction, n));// Ideal dielectric REFRACTION
    bool into = Vector3f::dot(nl, n) > 0;               // Ray from outside going in?
    double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = Vector3f::dot(nl, r.direction), cos2t;
    if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0)// Total internal reflection
        return obj->e + f * (radiance(reflRay, depth, Xi));
    Vector3f tdir = (r.direction * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).normalized();
    double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : Vector3f::dot(n,tdir));
    double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
    return obj->e + f * (depth > 2 ? (erand48(Xi) < P ?// Russian roulette
                                      radiance(reflRay, depth, Xi) * RP
                                                      : radiance(Ray(x, tdir), depth, Xi) * TP)
                                   : radiance(reflRay, depth, Xi) * Re + radiance(Ray(x, tdir), depth, Xi) * Tr);

}

int main(int argc, char *argv[]) {
    int w = 512, h = 384, samps = argc >= 2 ? atoi(argv[1]) / 4 : 5;// # samples

    double flength = 0;
    double aperture = 6;

    Ray cam(Vector3f(50, 52, 295.6), Vector3f(0, -0.042612, -1).normalized());      // cam pos, dir



    Vector3f cx = Vector3f(w * .5135 / h,0,0), cy =
            Vector3f::cross(cx, cam.direction).normalized() * .5135, r, *c = new Vector3f[w * h];

#pragma omp parallel for schedule(dynamic, 1) private(r)// OpenMP
    for (int y = 0; y < h; y++) {                       // Loop over image rows
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100. * y / (h - 1));
        for (unsigned short x = 0, Xi[3] = {0, 0, static_cast<unsigned short>(y * y * y)}; x < w; x++)// Loop cols
            for (int sy = 0, i = (h - y - 1) * w + x; sy < 2; sy++)                                 // 2x2 subpixel rows
                for (int sx = 0;
                     sx < 2; sx++, r = Vector3f()) {                                           // 2x2 subpixel cols
                    for (int s = 0; s < samps; s++) {

                        if (flength == 0) {
                            double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                            double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                            Vector3f d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                                         cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.direction;
                            r = r + radiance(Ray(cam.origin + d * 140, d.normalized()), 0, Xi) * (1. / samps);
                        } else {
                            double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                            double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                            Vector3f d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                                         cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.direction;
                            d = d.normalized() * flength;

                            Vector3f p;
                            do {
                                p = Vector3f(drand48(), drand48(), 0) * 2 - Vector3f(1, 1, 0);
                            } while (Vector3f::dot(p, p) >= 1);

                            Vector3f origin = cam.origin + p * aperture / 2;
                            Vector3f direction = d - p * aperture / 2;
                            r = r + radiance(Ray(origin, direction.normalized()), 0, Xi) * (1. / samps);
                        }


                    }// Camera rays are pushed ^^^^^ forward to start in interior
                    c[i] = c[i] + Vector3f(clamp(r.x()), clamp(r.y()), clamp(r.z())) * .25;
                }
    }
    FILE *f = fopen("news.ppm", "w");// Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i = 0; i < w * h; i++)
        fprintf(f, "%d %d %d ", toInt(c[i].x()), toInt(c[i].y()), toInt(c[i].z()));
    return 0;
}
