#pragma once

#include <vector>
#include <string>
#include "particles.h"

using namespace std;

void export_pointcloud(const vector<Particle> *particles, string path);
