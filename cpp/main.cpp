#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <QFile>
#include <QImage>
#include <QImageWriter>

#include "color.h"
#include "mymath.h"
#include "data.h"
#include "particles.h"
#include "image.h"
#include "export_ply.h"
#include "linking.h"

using namespace std;

//=============================================================================
// main
//=============================================================================

void connexion_cost(const float C,
                    const Particle *itp,
                    const vector<vector<vector<vector<vector<Particle*>>>>> &grid,
                    float *sum_da,
                    float *sum_proj_on_normal_h,
                    float *sum_w)
{
    *sum_da = 0.f;
    *sum_proj_on_normal_h = 0.f;
    *sum_w = 0.f;

    const float cos_n_h = cos(itp->h/D_HEADING * 2.f*M_PI + M_PI*0.5f);
    const float sin_n_h = sin(itp->h/D_HEADING * 2.f*M_PI + M_PI*0.5f);

    const int ix = (int)itp->x;
    const int iy = (int)itp->y;
    const int ih = (int)(itp->h - .5f);
    const int ia = (int)(itp->a - .5f);

    const float fh = (itp->h - .5f) - ih;
    const float fa = (itp->a - .5f) - ia;

    const int is = 4/(1<<DOWNSCALE);
    const float sigma = C * (float)is;
    const int sigma_crop = max(1, int(sigma*2.5f));

    for (int dh=0; dh<2; dh++)
    {
        const int ih2 = (ih+dh + D_HEADING) % D_HEADING;
        const float weight_h = dh==0 ? 1.f-fh : fh;

        for (int da=0; da<2; da++)
        {
            const int ia2 = max(0, min(D_ALTITUDE-1, ia+da));
            const float weight_a = da==0 ? 1.f-fa : fa;

            for (int dx=-sigma_crop; dx<=sigma_crop; dx++)
            {
                for (int dy=-sigma_crop; dy<=sigma_crop; dy++)
                {
                    // gaussian weight: TODO compute only once
                    const float gaussian_w = exp(-0.5f*(dx*dx + dy*dy)/(sigma*sigma));

                    if (0 <= ix+dx && ix+dx < D_WIDTH &&
                        0 <= iy+dy && iy+dy < D_HEIGHT)
                    {
                        for (auto itp2  = grid[ix+dx][iy+dy][ih2][ia2].begin();
                                  itp2 != grid[ix+dx][iy+dy][ih2][ia2].end(); itp2++)
                        {
                            float n2 = ((*itp2)->x - itp->x)*((*itp2)->x - itp->x) +
                                       ((*itp2)->y - itp->y)*((*itp2)->y - itp->y);
                            
                            if (n2 > 0.f)
                            {
                                const float w = weight_h * weight_a * gaussian_w * (*itp2)->w;

                                *sum_da               += w * ((*itp2)->a - itp->a);
                                *sum_proj_on_normal_h += w * (cos_n_h*((*itp2)->x - itp->x) +
                                                              sin_n_h*((*itp2)->y - itp->y)) / sqrt(n2);
                                *sum_w += w;
                            }
                        }
                    }
                }
            }
        }
    }
}

struct Params
{
    float A = 1.f;
    float B = 0.5f;
    float C = 2.f * 1.4f;
    float D = 1.f;
    float E = 1.f;
    float F = 1.f;

    string to_string()
    {
        return std::to_string(A) +"_"+
               std::to_string(B) +"_"+
               std::to_string(C) +"_"+
               std::to_string(D) +"_"+
               std::to_string(E) +"_"+
               std::to_string(F);
    }
};

void step0(vector<Particle> *particles, const Data &data, const float step_size)
{
    const float factors[3] = {16.f/(float)(1<<DOWNSCALE), 0.f, 1.f};

#pragma omp parallel for
    for (auto itp = particles->begin(); itp != particles->end(); itp++)
    {
        const float dx = factors[0] > 0.f ? data.get_derivative(itp->x, itp->y, itp->h, itp->a, 1) : 0.f;
        const float dy = factors[0] > 0.f ? data.get_derivative(itp->x, itp->y, itp->h, itp->a, 0) : 0.f;
        const float dh = factors[1] > 0.f ? data.get_derivative(itp->x, itp->y, itp->h, itp->a, 2) : 0.f;
        const float da = factors[2] > 0.f ? data.get_derivative(itp->x, itp->y, itp->h, itp->a, 3) : 0.f;

        const float cos_n_h = cos(itp->h/D_HEADING * 2.f*M_PI + M_PI*0.5f);
        const float sin_n_h = sin(itp->h/D_HEADING * 2.f*M_PI + M_PI*0.5f);
        float deriv_n = dx*cos_n_h + dy*sin_n_h;

        itp->x += factors[0] * (deriv_n >= 0.f ? step_size : -step_size) * cos_n_h;
        itp->y += factors[0] * (deriv_n >= 0.f ? step_size : -step_size) * sin_n_h;
        itp->h += factors[1] * (dh      >= 0.f ? step_size : -step_size);
        itp->a += factors[2] * (da      >= 0.f ? step_size : -step_size);

        while (itp->h < 0.f)        itp->h += D_HEADING;
        while (itp->h >= D_HEADING) itp->h -= D_HEADING;
        itp->a = clip(itp->a, 1e-6f, (float)D_ALTITUDE-1e-6f);
    }
}

void step1(vector<Particle> *particles0, vector<Particle> *particles, const Data &data, Params params, const float step_size)
{
    const float factors[3] = {16.f/(float)(1<<DOWNSCALE), 1.f, 2.f};

    //=====================================================================
    // compute grid
    //=====================================================================

    vector<vector<vector<vector<vector<Particle*>>>>> grid(D_WIDTH,
                                                           vector<vector<vector<vector<Particle*>>>>(D_HEIGHT,
                                                           vector<vector<vector<Particle*>>>(D_HEADING,
                                                           vector<vector<Particle*>>(D_ALTITUDE))));

    for (auto itp = particles->begin(); itp != particles->end(); itp++)
    {
        int ix = (int)itp->x;
        int iy = (int)itp->y;
        int ih = (int)itp->h;
        int ia = (int)itp->a;
        
        assert(0 <= ih && ih < D_HEADING);
        assert(0 <= ia && ia < D_ALTITUDE);

        if (0 <= ix && ix < D_WIDTH &&
            0 <= iy && iy < D_HEIGHT)
        {
            grid[ix][iy][ih][ia].push_back(&(*itp));
        }
    }

    //=====================================================================
    //
    //=====================================================================
    
    const size_t nb_particles = particles->size();

#pragma omp parallel for
    for (size_t p=0; p<nb_particles; p++)
    {
        Particle* itp  = &((*particles)[p]);
        Particle* itp0 = &((*particles0)[p]);

        const float dx = factors[0] > 0.f ? data.get_derivative(itp->x, itp->y, itp->h, itp->a, 1) : 0.f;
        const float dy = factors[0] > 0.f ? data.get_derivative(itp->x, itp->y, itp->h, itp->a, 0) : 0.f;
        const float dh = factors[1] > 0.f ? data.get_derivative(itp->x, itp->y, itp->h, itp->a, 2) : 0.f;
        const float da = factors[2] > 0.f ? data.get_derivative(itp->x, itp->y, itp->h, itp->a, 3) : 0.f;
        const float value = data.get_value(itp->x, itp->y, itp->h, itp->a);

        const float cos_n_h = cos(itp->h/D_HEADING * 2.f*M_PI + M_PI*0.5f);
        const float sin_n_h = sin(itp->h/D_HEADING * 2.f*M_PI + M_PI*0.5f);
        float deriv_n = dx*cos_n_h + dy*sin_n_h;

#if 0
        float d_dist_0_dx = itp0->x - itp->x;
        float d_dist_0_dy = itp0->y - itp->y;
        float d_dist_0_dn = d_dist_0_dx*cos_n_h + d_dist_0_dy*sin_n_h;
#else
        float d_dist_0_dn = 0.f;
#endif
        
        // connexion cost
        float continuity_a=0.f, continuity_xy=0.f, sum_w=0.f;
#if 0
        connexion_cost(params.C, itp, grid, &continuity_a, &continuity_xy, &sum_w);
#endif
        const float continuity_cst = 0.003f;

        const float deriv_n2 = deriv_n +
                               params.E * continuity_cst*continuity_xy +
                               params.F * value * d_dist_0_dn;
        itp->x += factors[0] * (deriv_n2 >= 0.f ? step_size : -step_size) * cos_n_h;
        itp->y += factors[0] * (deriv_n2 >= 0.f ? step_size : -step_size) * sin_n_h;
        
        
        itp->debug[0] = dh;
        itp->debug[1] = value*(itp->h0 - itp->h)*fabs(itp->h0 - itp->h);
#if 0
        const float dh_2 = dh + params.D*value*(itp->h0 - itp->h)*fabs(itp->h0 - itp->h);
        itp->h += factors[1] * (dh_2 >= 0.f ? step_size : -step_size);
#elif 0
        const float dh_2 = dh;
        itp->h += factors[1] * (dh_2 >= 0.f ? step_size : -step_size);
#endif

        itp->a += factors[2] * (da   >= 0.f ? step_size : -step_size);

        while (itp->h < 0.f)        itp->h += D_HEADING;
        while (itp->h >= D_HEADING) itp->h -= D_HEADING;

        itp->a = clip(itp->a, 1e-6f, (float)D_ALTITUDE-1e-6f);
    }
}



int main()
{
    const int grid_smoothing_iterations = 25; //45;
    const string particle_path =  "../results/"+ to_string(grid_smoothing_iterations) +".particles";

    //=========================================================================
    // 0) some input visualization
    //=========================================================================

#if 0
    {
        Data data("../data/data_11.f32", grid_smoothing_iterations);
        //Data data("../WebGL/data_11.f32", grid_smoothing_iterations);
        export_image(&data, "../results/vector_field.png");
    }
#endif

    //=========================================================================
    // 1) particle optimization
    //=========================================================================

    const float particle_threshold = 0.1f;

    if (!QFile::exists(particle_path.c_str()))
    {
        cout << particle_path << " NOT FOUND, computing..." << endl;
        const vector<float> param_vals = {0.f}; //{0.f, 0.1f, 1.f};

        for (size_t param_i=0; param_i<param_vals.size(); param_i++)
        {
            Params params;
            params.F = param_vals[param_i];

            // data
            Data data("../data/data_11.f32", grid_smoothing_iterations);

            // particles
            vector<Particle> particles0, particles;
            data.get_particles(&particles0, particle_threshold);
            data.get_particles(&particles,  particle_threshold);
            cout << "nb particles: " << particles.size() << endl;
            cout << D_HEIGHT << "x" << D_WIDTH << "x" << D_HEADING << "x" <<  D_ALTITUDE << endl;


            // optimize particles
            const int nb_iterations = 501;
            for (int it=0; it<nb_iterations; it++)
            {
                //const float step_size = 0.03f;
                const float step_size = 0.03f * pow(0.997f, (float)(it+1));
#if 0
                step0(&particles0, data, step_size);
#elif 1
                step1(&particles0, &particles, data, params, step_size);
#else
                step0(&particles0, data, step_size);
                step1(&particles0, &particles, data, params, step_size);
#endif
                if (it % 20 == 0 || it==nb_iterations-1)
                {
                    string filename = string("../results/")+ std::to_string(grid_smoothing_iterations) +"_"+ params.to_string();
#if 0
                    export_image(&data, &particles0, &particles0, filename +"_"+ to_string(it) +".png");
                    export_pointcloud(&particles0, filename +".ply");
#else
                    export_image(&data, &particles0, &particles, filename +"_"+ to_string(it) +".png");
                    //export_pointcloud(&particles, filename +".ply");
                    export_pointcloud(&particles, "../results/last.ply");
#endif
                }
                cout << it << endl;
            }

            Particle::save_particles(&particles, particle_path);
        }
    }

    //=========================================================================
    // 2) trajectory inference from particles
    //=========================================================================

#if 1
    {
        Data data("../data/data_11.f32", grid_smoothing_iterations);

        // particles before optimization
        vector<Particle> particles0;
        data.get_particles(&particles0, particle_threshold);

        // particles after optimization
        vector<Particle> particles;
        Particle::load_particles(&particles, particle_path);

        // skyways inference
        linking(&data, &particles, &particles0);
    }
#endif

    return 0;
}