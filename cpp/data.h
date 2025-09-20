#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "particles.h"


using namespace std;

//
#define DOWNSCALE     2

//
#define D0_HEIGHT  1024
#define D0_WIDTH   1024
#define D0_HEADING   12
#define D0_ALTITUDE  12

//
#define D_HEIGHT   (D0_HEIGHT/(1 << DOWNSCALE))
#define D_WIDTH    (D0_WIDTH /(1 << DOWNSCALE))
#define D_HEADING   D0_HEADING
#define D_ALTITUDE  D0_ALTITUDE

//
#if 0
#define CROP_Y0 ( 890/8)
#define CROP_Y1 (1755/8)
#define CROP_X0 ( 580/8)
#define CROP_X1 (1356/8)
#else
#define CROP_Y0 0
#define CROP_Y1 D_HEIGHT
#define CROP_X0 0
#define CROP_X1 D_WIDTH
#endif


class Data
{
    vector<float> data;
    
    // original
    const int f0y = D0_WIDTH * D0_HEADING * D0_ALTITUDE, fiy = D0_WIDTH * D0_HEADING * D_ALTITUDE, fy = D_WIDTH * D_HEADING * D_ALTITUDE;
    const int f0x =            D0_HEADING * D0_ALTITUDE, fix =            D0_HEADING * D_ALTITUDE, fx =           D_HEADING * D_ALTITUDE;
    const int f0h =                         D0_ALTITUDE, fih =                         D_ALTITUDE, fh =                       D_ALTITUDE;
    const int f0a =                                   1, fia =                                  1, fa =                                1;

public:
    Data(const string path, const int nb_smoothing_iterations);

    // create 1 particle per cell
    void get_particles(vector<Particle> *particles, const float threshold);

    float get(int x, int y, int h, int a) const;

    void prepare_interpolation(float x, float y, float p[4][4], float coords[2]) const;

    void prepare_interpolation(float x, float y, float a, float p[4][4][4], float coords[3]) const;

    void prepare_interpolation(float x, float y, float h, float a, float p[4][4][4][4], float coords[4]) const;

    float get_derivative(float x, float y, float h, float a, int dim, int order=1) const;
    
    float get_value(float x, float y, float h, float a) const;
};
