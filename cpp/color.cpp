#include <math.h>
#include "color.h"


RGB HSVtoRGB(const HSV& hsv)
{
    float h = hsv.h;
    float s = hsv.s;
    float v = hsv.v;

    RGB rgb;
    if (s == 0) {
        // If saturation is 0, the color is gray (R=G=B=V)
        rgb.r = rgb.g = rgb.b = v;
        return rgb;
    }

    float c = v * s;                  // Chroma
    float x = c * (1 - std::fabs(std::fmod(h / 60.0f, 2) - 1)); // Intermediate value
    float m = v - c;                  // Match value to bring to [0, 1]
    float r, g, b;

    if (h >= 0 && h < 60) {
        r = c; g = x; b = 0;
    } else if (h >= 60 && h < 120) {
        r = x; g = c; b = 0;
    } else if (h >= 120 && h < 180) {
        r = 0; g = c; b = x;
    } else if (h >= 180 && h < 240) {
        r = 0; g = x; b = c;
    } else if (h >= 240 && h < 300) {
        r = x; g = 0; b = c;
    } else { // h >= 300 && h < 360
        r = c; g = 0; b = x;
    }

    // Add the match value to shift to [0, 1]
    rgb.r = r + m;
    rgb.g = g + m;
    rgb.b = b + m;
    return rgb;
}

