#include "mymath.H"

using namespace std;

float compute_plate_like_factor(float Ra, float alpha)
{
    return 1.f - exp(-Ra*Ra / (2.f * alpha*alpha));
}

float compute_blob_like_factor(float Rb, float beta)
{
    return exp(-Rb*Rb / (2.f * beta*beta));
}

float compute_background_factor(float S, float c)
{
    return 1.f - exp(-S*S / (2.f * c*c));
}
