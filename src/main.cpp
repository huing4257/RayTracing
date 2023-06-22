#include <math.h>  // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdio.h> //        Remove "-fopenmp" for g++ version < 4.2
#include <stdlib.h>// Make : g++ -O3 -fopenmp smallpt.cpp -origin smallpt
#include "../include/vector.h"
#include "../include/sphere.hpp"
#include "../include/ray.hpp"
#include "../include/object3d.hpp"
#include "../include/hit.h"
#include "../include/plane.hpp"
#include "../include/triangle.hpp"
#include "../include/mesh.hpp"
#include "../include/revsurface.hpp"

BezierCurve pCurve(std::vector<Vec>{Vec(10, 0, 0),
                                    Vec(2, 1, 0),
                                    Vec(5, 10, 0),
                                    Vec(10, 20, 0)});

Mapping texture = Mapping("texture/wood.jpg");
NormalMapping normalTexture = NormalMapping("texture/Wall_n.png");

Object3D *objs[] = {
        new Plane(Vec(1, 0, 0), 0, Vec(0, 0, 0), Vec(.75, .75, .75), DIFF, &normalTexture),  //Left
        new Plane(Vec(-1, 0, 0), -100, Vec(0, 0, 0), Vec(.75, .75, .75), DIFF, &normalTexture),//Rght
        new Plane(Vec(0, 0, 1), 0, Vec(), Vec(.75, .75, .75), DIFF),        //Back
        new Plane(Vec(0, 0, -1), -170, Vec(), Vec(), DIFF, &normalTexture),              //Frnt
        new Plane(Vec(0, 1, 0), 0, Vec(), Vec(.75, .75, .75), DIFF, &texture),        //Botm
        new Plane(Vec(0, -1, 0), -81.6, Vec(), Vec(.75, .75, .75), DIFF),//Top
        new Sphere(600, Vec(50, 681.6 - .27, 81.6), Vec(12, 12, 12), Vec(), DIFF),   //Lite
        new Sphere(10,Vec(40,10,70),Vec(),Vec(.99,.99,.99),REFR)
        // new Mesh("mesh/bunny_200.obj",DIFF,Vec(0,0,0),Vec(.75,.75,.75))
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

Vec radiance(const Ray &r, int depth, unsigned short *Xi) {
    double t; // distance to intersection
    int id = 0;                            // id of intersected object
    Hit hit;
    if (!intersect(r, t, id, hit)) return {};// if miss, return black
    const Object3D *obj = objs[id];       // the hit object
    Vec n = hit.hit_normal, nl = n.dot(r.direction) < 0 ? n : n * -1;
    Vec x = hit.hit_pos, f = hit.hit_color;
    double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y
                                                        : f.z;// max refl
    if (++depth > 5) {
        if (erand48(Xi) < p) f = f * (1 / p);
        else
            return obj->e;  //R.R.
    }

    if (obj->refl == DIFF) {
        // Ideal DIFFUSE reflection
        double r1 = 2 * M_PI * erand48(Xi), r2 = erand48(Xi), r2s = sqrt(r2);
        Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm(), v = w % u;
        Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
        return obj->e + f.mult(radiance(Ray(x, d), depth, Xi));

    } else if (obj->refl == SPEC) {// Ideal SPECULAR reflection
        return obj->e + f.mult(radiance(Ray(x, r.direction - n * 2 * n.dot(r.direction)), depth, Xi));
    }

    Ray reflRay(x, r.direction - n * 2 * n.dot(r.direction));// Ideal dielectric REFRACTION
    bool into = n.dot(nl) > 0;               // Ray from outside going in?
    double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.direction.dot(nl), cos2t;
    if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0)// Total internal reflection
        return obj->e + f.mult(radiance(reflRay, depth, Xi));
    Vec tdir = (r.direction * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
    double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : tdir.dot(n));
    double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
    return obj->e + f.mult(depth > 2 ? (erand48(Xi) < P ?// Russian roulette
                                        radiance(reflRay, depth, Xi) * RP
                                                        : radiance(Ray(x, tdir), depth, Xi) * TP)
                                     : radiance(reflRay, depth, Xi) * Re + radiance(Ray(x, tdir), depth, Xi) * Tr);

}

int main(int argc, char *argv[]) {
    int w = 1024, h = 768, samps = argc == 2 ? atoi(argv[1]) / 4 : 50;// # samples

    Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm());      // cam pos, dir

    Vec cx = Vec(w * .5135 / h), cy = (cx % cam.direction).norm() * .5135, r, *c = new Vec[w * h];

#pragma omp parallel for schedule(dynamic, 1) private(r)// OpenMP
    for (int y = 0; y < h; y++) {                       // Loop over image rows
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100. * y / (h - 1));
        for (unsigned short x = 0, Xi[3] = {0, 0, static_cast<unsigned short>(y * y * y)}; x < w; x++)// Loop cols
            for (int sy = 0, i = (h - y - 1) * w + x; sy < 2; sy++)                                 // 2x2 subpixel rows
                for (int sx = 0;
                     sx < 2; sx++, r = Vec()) {                                           // 2x2 subpixel cols
                    for (int s = 0; s < samps; s++) {
                        double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                                cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.direction;
                        r = r + radiance(Ray(cam.origin + d * 140, d.norm()), 0, Xi) * (1. / samps);
                    }// Camera rays are pushed ^^^^^ forward to start in interior
                    c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
                }
    }
    FILE *f = fopen("news.ppm", "w");// Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i = 0; i < w * h; i++)
        fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
    system("python convert.py");
}
