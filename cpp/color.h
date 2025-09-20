#pragma once

struct HSV {
    float h; // Hue [0, 360]
    float s; // Saturation [0, 1]
    float v; // Value [0, 1]
};

struct RGB {
    float r; // Red [0, 1]
    float g; // Green [0, 1]
    float b; // Blue [0, 1]
};

RGB HSVtoRGB(const HSV& hsv);
