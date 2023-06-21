//
// Created by 胡 on 2023/6/21.
//

#ifndef RT_TEXTURE_H
#define RT_TEXTURE_H
#include "vector.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <iostream>

class Texture{
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
    Vec get_color(double u, double v) {
        int i = (int)(u * width);
        int j = (int)(v * height);
        if (i < 0) i = 0;
        if (j < 0) j = 0;
        if (i >= width) i = width - 1;
        if (j >= height) j = height - 1;
        int index = (i + j * width) * channel_num;
        double r = int(tex_data[index]) / 255.0;
        double g = int(tex_data[index+ 1]) / 255.0;
        double b = int(tex_data[index + 2]) / 255.0;
        return {r, g, b};
    }
    int width, height, channel_num;
    unsigned char *tex_data;

};
#endif //RT_TEXTURE_H
