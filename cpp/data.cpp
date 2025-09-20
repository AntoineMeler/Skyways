#include "data.h"
#include "mymath.h"


Data::Data(const string path, const int nb_smoothing_iterations)
{
    this->data.resize(D0_HEIGHT * D0_WIDTH * D0_HEADING * D0_ALTITUDE);

    ifstream input(path, ios::binary);
    input.read(reinterpret_cast<char*>(this->data.data()), D0_HEIGHT * D0_WIDTH * D0_HEADING * D0_ALTITUDE * sizeof(float));
    input.close();

    //=====================================================================
    // smooth grid along x & y axes
    //=====================================================================

    vector<float> data_line[2] = {vector<float>(max(D0_WIDTH, D0_HEIGHT) * D0_HEADING * D0_ALTITUDE),
                                  vector<float>(max(D0_WIDTH, D0_HEIGHT) * D0_HEADING * D0_ALTITUDE)};

    for (int direction=0; direction<2; direction++)
    {
        for (int i=0; i<(direction==0 ? D0_WIDTH : D0_HEIGHT); i++)
        {
            for (int j=0; j<(direction==0 ? D0_HEIGHT : D0_WIDTH); j++)
                for (int h=0; h<D0_HEADING; h++)
                    for (int a=0; a<D0_ALTITUDE; a++)
                        data_line[0][j*D0_HEADING*D0_ALTITUDE + h*D0_ALTITUDE + a] =
                            this->data[(direction == 0 ? j*this->f0y + i*this->f0x : i*this->f0y + j*this->f0x) + h*this->f0h + a*this->f0a];
            
            for (int smoothing_it=0; smoothing_it<nb_smoothing_iterations; smoothing_it++)
            {
                // smooth line/column
                for (int j=0; j<(direction==0 ? D0_HEIGHT : D0_WIDTH); j++)
                {
                    float *p_dst    =  &data_line[(smoothing_it+1)%2][j*D0_HEADING*D0_ALTITUDE];
                    float *p_src[3] = {&data_line[ smoothing_it%2   ][                                      max(0, j-1)*(D0_HEADING*D0_ALTITUDE)],
                                       &data_line[ smoothing_it%2   ][                                             j   *(D0_HEADING*D0_ALTITUDE)],
                                       &data_line[ smoothing_it%2   ][min((direction==0 ? D0_HEIGHT : D0_WIDTH)-1, j+1)*(D0_HEADING*D0_ALTITUDE)]};
                    
                    for (int ha=0; ha<D0_HEADING*D0_ALTITUDE; ha++)
                        *(p_dst++) = 0.25f * (*(p_src[0]++) + *(p_src[2]++)) + 0.5f * (*(p_src[1]++));
                }
            }

            for (int j=0; j<(direction==0 ? D0_HEIGHT : D0_WIDTH); j++)
                for (int h=0; h<D0_HEADING; h++)
                    for (int a=0; a<D0_ALTITUDE; a++)
                        this->data[(direction == 0 ? j*this->f0y + i*this->f0x : i*this->f0y + j*this->f0x) + h*this->f0h + a*this->f0a] =
                            data_line[nb_smoothing_iterations%2][j*D0_HEADING*D0_ALTITUDE + h*D0_ALTITUDE + a];
        }
    }
    
    //=====================================================================
    // downscale space (copy)
    //=====================================================================

    if (DOWNSCALE > 0)
    {
        vector<float> data_downscaled(D_HEIGHT * D_WIDTH * D_HEADING * D_ALTITUDE, 0.f);
        for (int y=0; y<D_HEIGHT; y++)
            for (int x=0; x<D_WIDTH; x++)
                for (int h=0; h<D_HEADING; h++)
                    for (int a=0; a<D_ALTITUDE; a++)
                        for (int dy=0; dy<1<<DOWNSCALE; dy++)
                            for (int dx=0; dx<1<<DOWNSCALE; dx++)
                                data_downscaled[y*this->fy + x*this->fx + h*this->fh + a*this->fa] += (data[(y*(1<<DOWNSCALE)+dy)*this->fiy + (x*(1<<DOWNSCALE)+dx)*this->fix + h*this->fih + a*this->fia]) / ((1<<DOWNSCALE)*(1<<DOWNSCALE));
        this->data = data_downscaled;
    }
}


void Data::get_particles(vector<Particle> *particles, const float threshold)
{
    vector<long int> grid(D_HEIGHT * D_WIDTH * D_HEADING * D_ALTITUDE, -1);

    size_t i_data = 0;
    for (int y=0; y<D_HEIGHT; y++)
        for (int x=0; x<D_WIDTH; x++)
            for (int h=0; h<D_HEADING; h++)
                for (int a=0; a<D_ALTITUDE; a++)
                {
                    if (this->data[i_data] > threshold)
                    {
                        if (CROP_X0 <= x && x < CROP_X1 && CROP_Y0 <= y && y < CROP_Y1)
                        {
                            const float w = this->data[i_data];
                            grid[i_data] = (long int)particles->size(); // id in particles[]
                            particles->push_back(Particle((float)x + 0.5f,
                                                          (float)y + 0.5f,
                                                          (float)h + 0.5f,
                                                          (float)a + 0.5f,
                                                          w));
                        }
                    }
                    i_data++;
                }
}

float Data::get(int x, int y, int h, int a) const
{
    return this->data[y*this->fy + x*this->fx + h*this->fh + a*this->fa];
}

void Data::prepare_interpolation(float x, float y, float p[4][4], float coords[2]) const
{
    int ix = (int)x;
    int iy = (int)y;

    coords[0] = y - iy;
    coords[1] = x - ix;

    int ipy[4], ipx[4];
    for (int dy=-1; dy<3; dy++) ipy[dy+1] = max(0, min(D_HEIGHT-1, iy+dy));
    for (int dx=-1; dx<3; dx++) ipx[dx+1] = max(0, min(D_WIDTH -1, ix+dx));

    for (int dy=0; dy<4; dy++)
    {
        for (int dx=0; dx<4; dx++)
        {
            p[dy][dx] = 0.f;
            for (int h=0; h<D_HEADING; h++)
                for (int a=0; a<D_ALTITUDE; a++)
                    p[dy][dx] += this->get(ipx[dx], ipy[dy], h, a);
        }
    }
}

void Data::prepare_interpolation(float x, float y, float a, float p[4][4][4], float coords[3]) const
{
    int ix = (int)x;
    int iy = (int)y;
    int ia = (int)a;

    coords[0] = y - iy;
    coords[1] = x - ix;
    coords[2] = a - ia;

    int ipy[4], ipx[4], ipa[4];
    for (int dy=-1; dy<3; dy++) ipy[dy+1] = max(0, min(D_HEIGHT  -1, iy+dy));
    for (int dx=-1; dx<3; dx++) ipx[dx+1] = max(0, min(D_WIDTH   -1, ix+dx));
    for (int da=-1; da<3; da++) ipa[da+1] = max(0, min(D_ALTITUDE-1, ia+da));

    for (int dy=0; dy<4; dy++)
    {
        for (int dx=0; dx<4; dx++)
        {
            for (int da=0; da<4; da++)
            {
                p[dy][dx][da] = 0.f;
                for (int h=0; h<D_HEADING; h++)
                    p[dy][dx][da] += this->get(ipx[dx], ipy[dy], h, ipa[da]);
            }
        }
    }
}

void Data::prepare_interpolation(float x, float y, float h, float a, float p[4][4][4][4], float coords[4]) const
{
    int ix = (int)x;
    int iy = (int)y;
    int ih = (int)h;
    int ia = (int)a;

    coords[0] = y - iy;
    coords[1] = x - ix;
    coords[2] = h - ih;
    coords[3] = a - ia;

    int ipy[4], ipx[4], iph[4], ipa[4];
    for (int dy=-1; dy<3; dy++) ipy[dy+1] = max(0, min(D_HEIGHT-1, iy+dy));
    for (int dx=-1; dx<3; dx++) ipx[dx+1] = max(0, min(D_WIDTH -1, ix+dx));
    for (int dh=-1; dh<3; dh++) iph[dh+1] = (ih+dh+D_HEADING)%D_HEADING;
    for (int da=-1; da<3; da++) ipa[da+1] = max(0, min(D_ALTITUDE-1, ia+da));

    for (int dy=0; dy<4; dy++)
        for (int dx=0; dx<4; dx++)
            for (int dh=0; dh<4; dh++)
                for (int da=0; da<4; da++)
                    p[dy][dx][dh][da] = this->get(ipx[dx], ipy[dy], iph[dh], ipa[da]);
}

float Data::get_derivative(float x, float y, float h, float a, int dim, int order) const
{
    // TODO: partager ça avec toutes les particules d'une même cellule
    float p[4][4][4][4];
    float coords[4];
    prepare_interpolation(x, y, h, a, p, coords);

    int derivative_order[4] = {dim==0 ? order : 0,
                               dim==1 ? order : 0,
                               dim==2 ? order : 0,
                               dim==3 ? order : 0};
    
    return nCubicInterpolate(4, (float*)p, coords, derivative_order);
}

float Data::get_value(float x, float y, float h, float a) const
{
    return get_derivative(x, y, h, a, 0, 0);
}