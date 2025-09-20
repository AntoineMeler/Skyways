#include <fstream>
#include <math.h>
#include "color.h"
#include "export_ply.h"
#include "data.h"


void export_pointcloud(const vector<Particle> *particles, string path)
{
    ofstream fout;
    fout.open(path);

    fout << "ply\n";
    fout << "format ascii 1.0\n";
    fout << "element vertex "<< particles->size() <<"\n";
    fout << "property float x\n";
    fout << "property float y\n";
    fout << "property float z\n";
    fout << "property uint8 red\n";
    fout << "property uint8 green\n";
    fout << "property uint8 blue\n";
    fout << "property uint8 alpha\n";
    fout << "end_header\n";

    for (auto itp=particles->begin(); itp != particles->end(); itp++)
    {
#if 1
        float angle = itp->h/D_HEADING * 360.f + 180.f;
        while (angle >= 360.f)
            angle -= 360.f;
        RGB rgb = HSVtoRGB({angle, 1., pow(itp->w, 0.5f)/2.f});
        int red   = max(0, min(255, int(256.f * rgb.r)));
        int green = max(0, min(255, int(256.f * rgb.g)));
        int blue  = max(0, min(255, int(256.f * rgb.b)));
        int alpha = 255; //max(0, min(255, int(1.f * itp->w)));
#else
        float v1 = itp->debug[0]/(itp->debug[0] + itp->debug[1]);
        float v2 = itp->debug[1]/(itp->debug[0] + itp->debug[1]);
        int red   = clip((int)(250.f*v2), 0, 255);
        int green = 150;
        int blue  = clip((int)(250.f*v1), 0, 255);
        int alpha = 255;
#endif

        fout << itp->x - D_WIDTH/2.f               << " " <<
                (D_HEIGHT - itp->y) - D_HEIGHT/2.f << " " <<
                4.f*itp->a                         << " " <<
                red << " " << green << " " << blue << " " << alpha << "\n";
    }
    fout.close();
}