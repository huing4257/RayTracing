#ifndef GROUP_H
#define GROUP_H


#include "object3d.hpp"
#include "ray.hpp"
#include "hit.hpp"
#include <iostream>
#include <vector>


// TODO: Implement Group - add data structure to store a list of Object*
class Group : public Object3D {

public:

    Group() {

    }

    explicit Group(int num_objects) {
        for (int i = 0; i < num_objects; ++i) {
            objects.push_back(nullptr);
        }
    }

    ~Group() override {

    }

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        bool isIntersect = false;
        for (auto object: objects) {
            if (object->intersect(r, h, tmin)) {
                isIntersect = true;
            }
        }
        return isIntersect;
    }

    void addObject(int index, Object3D *obj) {
        objects[index] = obj;
    }

    int getGroupSize() {
        return (int) objects.size();
    }

private:
    std::vector<Object3D *> objects;
};

#endif
	
