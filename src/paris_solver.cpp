#include "coords.hpp"
#include "isa.hpp"
#include "utils.hpp"
#include "solver.hpp"
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <cmath>
#include <iostream>

using namespace ballistics;

constexpr double g0 = 9.80665;

int ode_system(double t, const double y[], double f[], void *params) {
    BallisticParams *p = static_cast<BallisticParams*>(params);

    // y[0] — горизонтальная дистанция; y[1] — высота
    double x     = y[0];
    double height = y[1];
    double vx    = y[2];
    double vy    = y[3];

    // Плотность и температура из ISA
    double rho = ISA::get_density(height);
    double T   = ISA::get_temperature(height);

    // Число Маха
    constexpr double gamma_air = 1.4, R_air = 287.05;
    double a_sound = std::sqrt(gamma_air * R_air * T);
    double v = std::hypot(vx, vy);
    double M = v / a_sound;

    // Сила сопротивления
    double cD = drag_coefficient(M);
    double drag = 0.5 * rho * v*v * cD * p->area;

    // Ускорения
    double ax = - (drag / p->mass) * (vx / v);
    double ay = - g0               - (drag / p->mass) * (vy / v);

    f[0] = vx;
    f[1] = vy;
    f[2] = ax;
    f[3] = ay;
    return GSL_SUCCESS;
}

int main() {
    BallisticParams params{ 1210.0, 0.045 };

    // стартовые геодезические → Cartesian
    double lat0 = 49.9500, lon0 = 5.82, h0 = 335.0;
    GeodeticCoordinates geo{lon0, lat0, h0};
    CartesianCoordinates start = geod_to_cart(geo);

    // Заполняем вектор состояний
    double y[4];
    y[0] = start.x;  // горизонтальная
    y[1] = start.z;  // высота
    double speed0 = 1640.0, elev = 55.0 * M_PI/180.0;
    y[2] = speed0 * std::cos(elev);
    y[3] = speed0 * std::sin(elev);

    // GSL‑решатель
    gsl_odeiv2_system sys = { ode_system, nullptr, 4, &params };
    auto *d = gsl_odeiv2_driver_alloc_y_new(
        &sys, gsl_odeiv2_step_rkf45, 1e-3, 1e-6, 1e-6);

    double t = 0.0, dt = 0.001, t1 = 200.0;
    std::cout << "t(s)\tx(m)\ty(m)\tvx(m/s)\tvy(m/s)\n";
    for (int i = 1; i <= int(t1/dt); ++i) {
        double ti = i*dt;
        if (gsl_odeiv2_driver_apply(d, &t, ti, y) != GSL_SUCCESS) break;
        if (y[1] <= start.z) break;  // коснулись земли
        std::cout << t << '\t'
                  << y[0] << '\t' << y[1] << '\t'
                  << y[2] << '\t' << y[3] << "\n";
    }
    gsl_odeiv2_driver_free(d);
    return 0;
}
