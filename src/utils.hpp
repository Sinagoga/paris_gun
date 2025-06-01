#ifndef UTILS_HPP_
#define UTILS_HPP_

#include "coords.hpp"

namespace ballistics {

    struct NewtonsResult {
        double root;
        int iterations;
        bool converged;
    };

    double g(CartesianCoordinates cart_coords, double t);
    double dg(CartesianCoordinates cart_coords, double t);
    NewtonsResult newton_method(CartesianCoordinates cart_coords, double epsilon = 1e-10, int iters = 1000);
    double calculate_h(CartesianCoordinates cart_coords, double latitude); // latitude in radians

}

#endif