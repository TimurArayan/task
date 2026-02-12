#include <iostream>
#include <Eigen/Dense>
#include <cstring>

#define REAL double
#define VOID void
extern "C" {
#include "triangle.h"
}

int main() {
    Eigen::Vector2d v(10, 20);
    std::cout << "Eigen: " << (v * 0.5).transpose() << std::endl;

   
    struct triangulateio in, out;
    std::memset(&in, 0, sizeof(in));
    std::memset(&out, 0, sizeof(out));

    REAL c[] = {0,0, 1,0, 1,1, 0,1};
    in.numberofpoints = 4;
    in.pointlist = c;

    triangulate((char*)"pqQ", &in, &out, nullptr);

    std::cout << "Triangle nodes: " << out.numberofpoints << std::endl;

    free(out.pointlist);
    free(out.trianglelist);
    return 0;
}