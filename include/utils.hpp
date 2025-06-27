#ifndef UTILS_HPP_
#define UTILS_HPP_

#include "coords.hpp"
#include <vector>

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
    double dmse(std::vector<double> M, std::vector<double> c_d_true, std::vector<double> c_d_pred);
    double get_a(std::vector<double> M, std::vector<double> C_d, double eta = 0.01, int iter = 1000);

}

#endif