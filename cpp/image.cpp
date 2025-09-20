
#include "mymath.h"

#include <QImage>
#include <QImageWriter>

#include "color.h"
#include "image.h"
#include "mymath.h"

#define RENDERING_UPSCALE 3


float divide(float a, float b)
{
    return b==0.f ? 0.f : a/b;
}


void export_image(const Data *data, string path)
{
#if 1
    vector<vector<vector<float>>> img_data(D_WIDTH, vector<vector<float>>(D_HEIGHT, vector<float>(3, 0.f)));

    for (int x=0; x<D_WIDTH; x++)
    {
        for (int y=0; y<D_HEIGHT; y++)
        {
            for (int a=0; a<D_ALTITUDE; a++)
            {
                //=============================================================
                // hessian
                //=============================================================

                float p[4][4][4];
                float coords[3];
                data->prepare_interpolation((float)x + .5f,
                                            (float)y + .5f,
                                            (float)a + .5f, p, coords);

                int derivative_order_dxdx[3] = {0, 2, 0};
                int derivative_order_dydy[3] = {2, 0, 0};
                int derivative_order_dxdy[3] = {1, 1, 0};
                int derivative_order_dxdz[3] = {0, 1, 1};
                int derivative_order_dydz[3] = {1, 0, 1};
                int derivative_order_dzdz[3] = {0, 0, 2};

                float dxdx = nCubicInterpolate(3, (float*)p, coords, derivative_order_dxdx);
                float dydy = nCubicInterpolate(3, (float*)p, coords, derivative_order_dydy);
                float dxdy = nCubicInterpolate(3, (float*)p, coords, derivative_order_dxdy);
                float dxdz = nCubicInterpolate(3, (float*)p, coords, derivative_order_dxdz);
                float dydz = nCubicInterpolate(3, (float*)p, coords, derivative_order_dydz);
                float dzdz = nCubicInterpolate(3, (float*)p, coords, derivative_order_dzdz);

                arma::mat H = {{dxdx, dxdy, dxdz},
                               {dxdy, dydy, dydz},
                               {dxdz, dydz, dzdz}};
                
                arma::vec eigval;
                arma::mat eigvec;
                arma::eig_sym(eigval, eigvec, H);

                //=============================================================
                // sort eigval by absolute value
                //=============================================================

                std::vector<std::pair<float, int>> indexed_eigval;
                for (size_t i=0; i<eigval.size(); i++)
                    indexed_eigval.emplace_back(eigval[i], i);
                
                std::sort(indexed_eigval.begin(), indexed_eigval.end(), [](const std::pair<float, int>& a, const std::pair<float, int>& b) { return std::abs(a.first) < std::abs(b.first); });

                std::vector<arma::vec> sorted_eigvec;
                std::vector<float> sorted_eigval;
                for (const auto& pair : indexed_eigval) {
                    sorted_eigval.push_back(pair.first);
                    sorted_eigvec.push_back(eigvec.col(pair.second));
                }
                //=============================================================
                // frangi: https://github.com/ellisdg/frangi3d/blob/master/frangi/frangi.py#L51
                //=============================================================

                float Ra = divide(fabs(sorted_eigval[1]), fabs(sorted_eigval[2])); // doit être grand (proche de 1)
                float Rb = divide(fabs(sorted_eigval[0]), sqrt(fabs(sorted_eigval[1] * sorted_eigval[2]))); // doit être petit
                float S  = sqrt(sorted_eigval[0]*sorted_eigval[0] +
                                sorted_eigval[1]*sorted_eigval[1] +
                                sorted_eigval[2]*sorted_eigval[2]);
                float alpha = .002f, beta = 0.5f, frangi_c = .5f;
                float plate = compute_plate_like_factor(Ra, alpha);
                float blob  = compute_blob_like_factor(Rb, beta);
                float background = compute_background_factor(S, frangi_c);
                float voxel_data = plate * blob * background;
                
                if (sorted_eigval[1] > 0.f || sorted_eigval[2] > 0.f)
                    voxel_data = 0.f;
#if 0
                img_data[x][y][0] += 0.05*plate;
                img_data[x][y][1] += 0.05*blob;
                img_data[x][y][2] += 0.05*background;
#else
                float v = 0.5 * voxel_data;
                img_data[x][y][0] += v;
                img_data[x][y][1] += v * sorted_eigvec[2][0];
                img_data[x][y][2] += v * sorted_eigvec[2][1];
#endif
            }
        }
    }

    QImage img(D_WIDTH, D_HEIGHT, QImage::Format_ARGB32);
    img.fill(QColor(0, 0, 0, 255));

    for (int ix=0; ix<D_WIDTH; ix++)
    {
        for (int iy=0; iy<D_HEIGHT; iy++)
        {
#if 0
            int red   = max(0, min(255, int(256.f * img_data[ix][iy][0])));
            int green = max(0, min(255, int(256.f * img_data[ix][iy][1])));
            int blue  = max(0, min(255, int(256.f * img_data[ix][iy][2])));
#else
            float angle = atan2(img_data[ix][iy][2], img_data[ix][iy][1])/M_PI * 180.;
            while (angle >= 360.f) angle -= 360.f;
            while (angle <    0.f) angle += 360.f;
            
            RGB rgb = HSVtoRGB({angle, 1., img_data[ix][iy][0]});
           
            int red   = max(0, min(255, int(256.f * rgb.r)));
            int green = max(0, min(255, int(256.f * rgb.g)));
            int blue  = max(0, min(255, int(256.f * rgb.b)));
#endif
            QColor color = QColor(red, green, blue, 255);
            
            img.setPixelColor(ix, iy, color);
        }
    }
    
    QImageWriter writer(QString(path.c_str()));
    writer.write(img);
#else
    const float step = 0.5f;
    vector<vector<float>> img_data(D_WIDTH, vector<float>(D_HEIGHT, 0.f));
    for (int p=0; p<100000; p++)
    {
        float x = (float)(rand()%D_WIDTH);
        float y = (float)(rand()%D_HEIGHT);
        float color = (float)(rand()%256)/255.f;

        for (int l=0; l<200; l++)
        {
            float p[4][4];
            float coords[2];
            data->prepare_interpolation(x, y, p, coords);

            int derivative_order_dx[2] = {0, 1};
            int derivative_order_dy[2] = {1, 0};

            float dx = nCubicInterpolate(2, (float*)p, coords, derivative_order_dx);
            float dy = nCubicInterpolate(2, (float*)p, coords, derivative_order_dy);
            float n2 = dx*dx + dy*dy;

            if (n2 > 0.f)
            {
                float rdx = -dy;
                float rdy = dx;
                x += step * rdx/sqrt(n2);
                y += step * rdy/sqrt(n2);
            }

            img_data[clip((int)x, 0, D_WIDTH-1)][clip((int)y, 0, D_HEIGHT-1)] = color;
        }
    }
    

    QImage img(D_WIDTH, D_HEIGHT, QImage::Format_ARGB32);
    img.fill(QColor(0, 0, 0, 255));

    for (int ix=0; ix<D_WIDTH; ix++)
    {
        for (int iy=0; iy<D_HEIGHT; iy++)
        {
            int vi = clip((int)(256.f * img_data[ix][iy]), 0, 255);
            QColor color = QColor(vi, vi, vi, 255);
            img.setPixelColor(ix, iy, color);
        }
    }
    
    QImageWriter writer(QString(path.c_str()));
    writer.write(img);
#endif
}


void save_image(const vector<vector<vector<float>>> &img_data, string path)
{
    QImage img(D_WIDTH*(1<<RENDERING_UPSCALE), D_HEIGHT*(1<<RENDERING_UPSCALE), QImage::Format_ARGB32);
    img.fill(QColor(0, 0, 0, 255));

    for (int ix=0; ix<D_WIDTH*(1<<RENDERING_UPSCALE); ix++)
    {
        for (int iy=0; iy<D_HEIGHT*(1<<RENDERING_UPSCALE); iy++)
        {
            if (img_data[ix][iy][0] > 0.f)
            {
#if 1
                float angle = (atan2(img_data[ix][iy][2], img_data[ix][iy][1]) + M_PI) * 180./M_PI;
                RGB rgb = HSVtoRGB({angle, 1., min(2.f * 0.075f*pow(img_data[ix][iy][0], 0.4f), 1.f)});
                QColor color = QColor(max(0, min(255, int(256.f * rgb.r))),
                                      max(0, min(255, int(256.f * rgb.g))),
                                      max(0, min(255, int(256.f * rgb.b))), 255);
#elif 1
                float v1 = img_data[ix][iy][3]/(img_data[ix][iy][3] + img_data[ix][iy][4]);
                float v2 = img_data[ix][iy][4]/(img_data[ix][iy][3] + img_data[ix][iy][4]);
                QColor color = QColor(clip((int)(250.f*v2), 0, 255),
                                      150,
                                      clip((int)(250.f*v1), 0, 255),
                                      255);
#else
                float v = 3.f * img_data[ix][iy][5];
                int vi = clip((int)(10.f*v), 0, 255);
                QColor color = QColor(vi, vi, vi, 255);
#endif
                img.setPixelColor(ix, iy, color);
            }
        }
    }
    
    QImageWriter writer(QString(path.c_str()));
    writer.write(img);
}


void export_image(const Data *data, vector<Particle> *particles0, vector<Particle> *particles, string path)
{
    vector<vector<vector<float>>> img_data(D_WIDTH*(1<<RENDERING_UPSCALE), vector<vector<float>>(D_HEIGHT*(1<<RENDERING_UPSCALE), vector<float>(5, 0.0f)));

    for (auto itp0=particles0->begin(), itp=particles->begin(); itp != particles->end(); itp0++, itp++)
    {
#if 1
        float x = (itp->x - 0.5)*(float)(1<<RENDERING_UPSCALE);
        float y = (itp->y - 0.5)*(float)(1<<RENDERING_UPSCALE);
#else
        float x = (itp0->x - 0.5)*(float)(1<<RENDERING_UPSCALE);
        float y = (itp0->y - 0.5)*(float)(1<<RENDERING_UPSCALE);
#endif
        int ix = int(x);
        int iy = int(y);
        float fx = x - ix;
        float fy = y - iy;
        float interp[2][2] = {{(1.f-fx)*(1.f-fy), (1.f-fx)*fy},
                              {     fx *(1.f-fy),      fx *fy}};

        for (int dx=0; dx<2; dx++)
            for (int dy=0; dy<2; dy++)
            {
                int dst_x = clip(ix+dx, 0, D_WIDTH  * (1<<RENDERING_UPSCALE)-1);
                int dst_y = clip(iy+dy, 0, D_HEIGHT * (1<<RENDERING_UPSCALE)-1);
                img_data[dst_x][dst_y][0] += interp[dx][dy] * itp->w;
                img_data[dst_x][dst_y][1] += interp[dx][dy] * itp->w * cos(itp->h/D_HEADING * 2.f*M_PI);
                img_data[dst_x][dst_y][2] += interp[dx][dy] * itp->w * sin(itp->h/D_HEADING * 2.f*M_PI);
                img_data[dst_x][dst_y][3] += interp[dx][dy] * fabs(itp->debug[0]);
                img_data[dst_x][dst_y][4] += interp[dx][dy] * fabs(itp->debug[1]);
            }
    }

    save_image(img_data, path);
}