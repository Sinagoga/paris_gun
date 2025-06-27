#include "adr.hpp"
#define _USE_MATH_DEFINES
#include <cmath>

namespace ballistics {

    double c_d(double m, double a, double b = 0.5) {
        return a / pow(m*m - 1, b);
    }

    double s(double diameter) {
        return (M_PI*diameter*diameter)/4; 
    }

    double R(double dens, double v, double s, double c_d) {
        return 1/2 * dens * v*v * s * c_d;
    }

}