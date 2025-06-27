#include "coords.hpp"
#include "utils.hpp"
#include "adr.hpp"
#include <cmath>
#include <iostream>

namespace ballistics {

    double g(CartesianCoordinates cart_coords, double t) {
        double u = sqrt(cart_coords.x*cart_coords.x + cart_coords.y*cart_coords.y);
        double z_a = cart_coords.z/R_eq;
        double u_a = u/R_eq;
        double c = (1-F)*(1-F);

        double g_t = z_a / (u_a - (F*(2-F))/sqrt(1 + c*t*t));

        return g_t;
    }

    double dg(CartesianCoordinates cart_coords, double t) {
        double u = sqrt(cart_coords.x*cart_coords.x + cart_coords.y*cart_coords.y);
        double u_a = u/R_eq;
        double z_a = cart_coords.z/R_eq;
        double c = (1-F)*(1-F);

        double psi = u_a - (F*(2 - F))/sqrt(1+c*t*t);
        double dpsi = (F*(2 - F)) * c*t/pow(1+c*t*t, 3/2);

        double dg_t = -(z_a*dpsi)/(psi*psi);

        return dg_t;
    }

    NewtonsResult newton_method(CartesianCoordinates cart_coords, double epsilon, int iters) {
        double u = sqrt(cart_coords.x*cart_coords.x + cart_coords.y*cart_coords.y);
        double t = cart_coords.z/u;

        for (int i = 0; i < iters; i++) {
            double t_new = t - (g(cart_coords, t) - t)/(dg(cart_coords, t) - 1);
            if (fabs(t_new - t) < epsilon) {
                return {t_new, i + 1, true};
            }
            t = t_new;
        }

        return {t, iters, false};
    }

    double calculate_h(CartesianCoordinates cart_coords, double latitude) {
        double u = sqrt(cart_coords.x*cart_coords.x + cart_coords.y*cart_coords.y);

        double h = (1/cos(latitude)) * (u - (R_eq*R_eq)/sqrt(R_eq*R_eq + R_pol*R_pol*tan(latitude)*tan(latitude)));

        return h;
    }

    
    double dmse(std::vector<double> M, std::vector<double> c_d_true, std::vector<double> c_d_pred) {
        int n = std::size(c_d_true);
        double sq_sum = 0;

        for (int i = 0; i < n; i++) {
            sq_sum += -2*(c_d_true[i]-c_d_pred[i])*(1/sqrt(M[i]*M[i]-1));
        }

        return sq_sum/n;
    }

    double get_a(std::vector<double> M, std::vector<double> C_d, double eta = 0.01, int iter = 1000) {
        double a = 1;

        for (int i = 0; i < iter; i++) {
            std::vector<double> c_d_pred = {};
            for (int j = 0; j < std::size(M); j++) {
                double pred = c_d(M[i], a, 0.5);
                c_d_pred.push_back(pred);
            }
            double dl = dmse(M, C_d, c_d_pred);
            a = a - eta * dl;
        }

        return a;
    }
}