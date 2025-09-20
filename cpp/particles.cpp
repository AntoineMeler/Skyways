#include <math.h>

#include "particles.h"


Particle::Particle() {}

Particle::Particle(float x, float y, float heading, float altitude, float weight)
{
    this->x = x;
    this->y = y;
    this->h = heading;
    this->a = altitude;
    this->w = weight;
    this->h0 = heading;
}

arma::vec Particle::get_xyah() const
{
    arma::vec pos = {this->x,
                     this->y,
                     this->a,
                     this->h};
    return pos;
}

void Particle::save(std::ofstream &out)
{
    out.write(reinterpret_cast<char*>(&x),  sizeof(float));
    out.write(reinterpret_cast<char*>(&y),  sizeof(float));
    out.write(reinterpret_cast<char*>(&h),  sizeof(float));
    out.write(reinterpret_cast<char*>(&a),  sizeof(float));
    out.write(reinterpret_cast<char*>(&w),  sizeof(float));
    out.write(reinterpret_cast<char*>(&h0), sizeof(float));
}

void Particle::load(std::ifstream &input)
{
    input.read(reinterpret_cast<char*>(&x),  sizeof(float));
    input.read(reinterpret_cast<char*>(&y),  sizeof(float));
    input.read(reinterpret_cast<char*>(&h),  sizeof(float));
    input.read(reinterpret_cast<char*>(&a),  sizeof(float));
    input.read(reinterpret_cast<char*>(&w),  sizeof(float));
    input.read(reinterpret_cast<char*>(&h0), sizeof(float));
}


void Particle::save_particles(vector<Particle> *particles, string path)
{
    size_t nb = particles->size();

    std::ofstream out(path, std::ios::out | std::ios::binary);

    out.write(reinterpret_cast<char*>(&nb),  sizeof(size_t));
    for (auto itp=particles->begin(); itp != particles->end(); itp++)
        itp->save(out);

    out.close();
}

void Particle::load_particles(vector<Particle> *particles, string path)
{
    size_t nb;

    std::ifstream input(path, std::ios::in | std::ios::binary);

    input.read(reinterpret_cast<char*>(&nb), sizeof(size_t));
    particles->resize(nb);
    for (auto itp=particles->begin(); itp != particles->end(); itp++)
        itp->load(input);

    input.close();
}