#include <iostream>
#include <cstdlib>

extern "C" {
    #define VOID void
    #define REAL double
    #include "triangle.h"
}

int main() {
    struct triangulateio in, out;

    in.pointlist = (REAL *) malloc(3 * 2 * sizeof(REAL));
    in.numberofpoints = 3;
    
    in.pointlist[0] = 0; in.pointlist[1] = 0;
    in.pointlist[2] = 1; in.pointlist[3] = 0;
    in.pointlist[4] = 0; in.pointlist[5] = 1;

    in.numberofpointattributes = 0;
    in.pointmarkerlist = nullptr;
    in.numberofsegments = 0;
    in.numberofholes = 0;
    in.numberofregions = 0;

    out.pointlist = nullptr;
    out.trianglelist = nullptr;

    triangulate((char*)"zQ", &in, &out, nullptr);

    std::cout << "Triangle OK! Nodes: " << out.numberofpoints 
              << ", Triangles: " << out.numberoftriangles << std::endl;

    free(in.pointlist);
    return 0;
}