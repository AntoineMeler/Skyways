#pragma once

#include <vector>
#include "data.h"

using namespace std;

void export_image(const Data *data, string path);
void save_image(const vector<vector<vector<float>>> &img_data, string path);
void export_image(const Data *data, vector<Particle> *particles0, vector<Particle> *particles, string path);
