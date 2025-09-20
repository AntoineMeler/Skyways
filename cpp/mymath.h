#pragma once
#include <cassert>
#include <armadillo>

using namespace std;

// https://www.paulinternet.nl/?page=bicubic


template<typename T>
T clip(T i, T minv, T maxv)
{
    return max(minv, min(maxv, i));
}

float compute_plate_like_factor(float Ra, float alpha);
float compute_blob_like_factor(float Rb, float beta);
float compute_background_factor(float S, float c);

template<typename T>
T cubicInterpolate(const T p[4], const T x, const int derivative_order)
{
    if (derivative_order == 0)
    {
        return p[1] + T(0.5) * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
    }
    else
    {
        const T a = 3.0*(p[1] - p[2]) + p[3] - p[0];
        const T b = 2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3];

        if (derivative_order == 1)
        {
            const T c = p[2] - p[0];
            return 1.5*a*x*x + b*x + 0.5*c; // 0.5*(3*a*x*x + 2*b*x + c)
        }
        else
        {
            return 3.*a*x + b;
        }
    }
}

template<typename T>
T nCubicInterpolate (int n, T* p, const T coordinates[], const int derivative_order[])
{
	assert(n > 0);
	if (n == 1) {
        return cubicInterpolate(p, *coordinates, *derivative_order);
	} else {
		T arr[4];
		int skip = 1 << (n - 1) * 2;
		arr[0] = nCubicInterpolate(n-1, p+0*skip, coordinates+1, derivative_order+1);
		arr[1] = nCubicInterpolate(n-1, p+1*skip, coordinates+1, derivative_order+1);
		arr[2] = nCubicInterpolate(n-1, p+2*skip, coordinates+1, derivative_order+1);
		arr[3] = nCubicInterpolate(n-1, p+3*skip, coordinates+1, derivative_order+1);
        return cubicInterpolate(arr, *coordinates, *derivative_order);
	}
}
