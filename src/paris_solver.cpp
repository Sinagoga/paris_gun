#include "coords.hpp"
#include "isa.hpp"
#include "utils.hpp"
#include "solver.hpp"
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <cmath>
#include <iostream>
#include "adr.hpp"

using namespace ballistics;

constexpr double g0 = 9.80665;

int ode_system(double t, const double y[], double f[], void *params) {
    BallisticParams *p = static_cast<BallisticParams *>(params);

    // y[0] — горизонтальная дистанция; y[1] — высота
    double x = y[0];
    double height = y[1];
    double vx = y[2];
    double vy = y[3];

    // Плотность и температура из ISA
    double rho = ISA::get_density(height);
    double T = ISA::get_temperature(height);

    // Число Маха
    constexpr double gamma_air = 1.4, R_air = 287.05;
    double a_sound = std::sqrt(gamma_air * R_air * T);
    double v = std::hypot(vx, vy);
    double M = v / a_sound;

    // Сила сопротивления
//    double cD = drag_coefficient(M);
    std::vector<double> M_ki = {1.1, 1.4, 2.2, 3, 4.2};
    std::vector<double> C_ki = {0.22, 0.2, 0.14, 0.1, 0.06};
//    a = ballistics::get_a(M_ki,C_ki)
    double cD = ballistics::c_d(M, 0.0856307,0.6);

    double drag = 0.5 * rho * v * v * cD * p->area;

    // Ускорения
    double ax = -(drag / p->mass) * (vx / v);
    double ay = -g0 - (drag / p->mass) * (vy / v);

    f[0] = vx;
    f[1] = vy;
    f[2] = ax;
    f[3] = ay;
    return GSL_SUCCESS;
}

int main() {



    double start_angle = 30.0;
    int res_angle = 0;
    double res_x = 0;
    for (int j = 40; j <= 75; j ++){
        start_angle += 1.0;
        // Заполняем вектор состояний
        BallisticParams params{1210.0, 0.045};

        // стартовые геодезические → Cartesian
        double lat0 = 49.9500, lon0 = 5.82, h0 = 335.0;
        GeodeticCoordinates geo{lon0, lat0, h0};
        CartesianCoordinates start = geod_to_cart(geo);
        double y[4];
        y[0] = start.x;  // горизонтальная
        y[1] = start.z;  // высота
        double speed0 = 1640.0, elev = start_angle * M_PI / 180.0;
        y[2] = speed0 * std::cos(elev);
        y[3] = speed0 * std::sin(elev);

        // GSL‑решатель
        gsl_odeiv2_system sys = {ode_system, nullptr, 4, &params};
        auto *d = gsl_odeiv2_driver_alloc_y_new(
                &sys, gsl_odeiv2_step_rkf45, 1e-3, 1e-6, 1e-6);
        double t = 0.0, dt = 0.1, t1 = 2000.0;
        //std::cout << "t(s)\tx(m)\ty(m)\tvx(m/s)\tvy(m/s)\n";
        double s = y[0] + 12e4;
        for (int i = 1; i <= int(t1 / dt); ++i) {
            double ti = i * dt;
            if (gsl_odeiv2_driver_apply(d, &t, ti, y) != GSL_SUCCESS) break;
            if (y[1] <= start.z) break;  // коснулись земли
        }
        double delta = (y[0]-s)/1000;
        std::cout<<"Угол возывшения (градусы) "<<j<<" | Дальность выстрела "<<delta<<"км\n";
        if (delta > res_x){
            res_x = delta;
            res_angle = j;
        }
        gsl_odeiv2_driver_free(d);
    }
    std::cout<<"\n\n\n"<<"Оптимальный угол возвышения - "<<res_angle<<"| Дальность выстрела - "<<res_x<<"км\n";
    return 0;
}
