#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <armadillo>

using namespace std;


class Particle
{
public:
    float x, y, h, a, w;
    float h0;
    int line_id = -1;
    float debug[2];

    Particle();
    Particle(float x, float y, float heading, float altitude, float weight);

    arma::vec get_xyah() const;

    void save(std::ofstream &out);
    void load(std::ifstream &input);

    static void save_particles(vector<Particle> *particles, string path);
    static void load_particles(vector<Particle> *particles, string path);
};