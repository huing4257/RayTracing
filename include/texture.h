//
// Created by 胡 on 2023/6/21.
//

#ifndef RT_TEXTURE_H
#define RT_TEXTURE_H

#include "vecmath.h"

#define STB_IMAGE_IMPLEMENTATION

#include "stb_image.h"
#include <iostream>

class Texture {
public:
    Texture(const char *filename) {
        unsigned char *data = stbi_load(filename, &width, &height, &channel_num, 0);
        if (data == nullptr) {
            std::cerr << "ERROR: could not load texture image file \"" << filename << "\"." << std::endl;
            exit(1);
        }
        tex_data = data;
    }

    ~Texture() = default;

    virtual void change_hit(double u, double v,Hit &hit) = 0;

    int width, height, channel_num;
    unsigned char *tex_data;

};

class Mapping : public Texture {
public:
    Mapping(const char *filename) : Texture(filename) {};
    void change_hit(double u, double v,Hit &hit) override {
        int i = (int) (u * width);
        int j = (int) (v * height);
        if (i < 0) i = 0;
        if (j < 0) j = 0;
        if (i >= width) i = width - 1;
        if (j >= height) j = height - 1;
        int index = (i + j * width) * channel_num;
        double r = int(tex_data[index]) / 255.0;
        double g = int(tex_data[index + 1]) / 255.0;
        double b = int(tex_data[index + 2]) / 255.0;
        hit.hit_color = Vector3f(r, g, b);
    }
};

class NormalMapping : public Texture {
public:
    NormalMapping(const char *filename) : Texture(filename) {};
    void change_hit(double u, double v,Hit &hit) override {
        int i = (int) (u * width);
        int j = (int) (v * height);
        if (i < 0) i = 0;
        if (j < 0) j = 0;
        if (i >= width) i = width - 1;
        if (j >= height) j = height - 1;
        int index = (i + j * width) * channel_num;
        // change the normal vector of hit
        double r = int(tex_data[index]) / 255.0;
        double g = int(tex_data[index + 1]) / 255.0;
        double b = int(tex_data[index + 2]) / 255.0;
        hit.hit_normal = hit.hit_normal + Vector3f(r - 0.5, g - 0.5, b - 0.5);
    }
};

#endif //RT_TEXTURE_H
