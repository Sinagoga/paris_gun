#include <iostream>
#include "utils.hpp"
#include "coords.hpp"
#include "isa.hpp"
// #include "isa_temperature.hpp" 
// #include "isa_pressure.hpp" 
#define _USE_MATH_DEFINES
#include <cmath>

using namespace ballistics;
int main() {
    const double test_heights[] = {0.0, 5000.0, 11000.0, 15000.0, 20000.0};

    const double etalon_T[]      = {288.15, 255.65, 216.65, 216.65, 216.65};

    std::cout << "ISA temperature checks:\n";
    for (int i = 0; i < 5; ++i) {
        double h = test_heights[i];
        double T = ISA::get_temperature(h);
        std::cout << "  h=" << h/1000 << " km: "
                  << T << " K (expected ~" << etalon_T[i] << ")"
                  << (std::fabs(T - etalon_T[i]) < 1e-6 ? " +" : " -")
                  << "\n";
    }
    std::cout << "\n";

    const double etalon_P[]      = {101325.0, 54019.0, 22632.0, 12041.0, 5474.0};

    std::cout << "ISA pressure checks:\n";
    for (int i = 0; i < 5; ++i) {
        double h = test_heights[i];
        double p = ISA::get_pressure(h);
        std::cout << "  h=" << h/1000 << " km: "
                << p << " Pa (expected ~" << etalon_P[i] << ")"
                << (std::fabs(p - etalon_P[i]) / etalon_P[i] < 0.01 ? " +" : " -")
                << "\n";
    }
    std::cout << "\n";

    const double etalon_D[]      = {1.225, 0.736, 0.364, 0.194, 0.088};

    std::cout << "ISA density checks:\n";
    for (int i = 0; i < 5; ++i) {
        double h = test_heights[i];
        double D = ISA::get_density(h);
        std::cout << "  h=" << h/1000 << " km: "
                  << D << " K (expected ~" << etalon_D[i] << ")"
                  << (std::fabs(D - etalon_D[i]) < 1e-3 ? " +" : " -")
                  << "\n";
    }
    std::cout << "\n";
    
    GeodeticCoordinates moscow;
    moscow.latitude = 55*M_PI / 180.0;
    moscow.longitude = 37*M_PI / 180.0;
    moscow.altitude = 124.0;

    std::cout << "Moscow geod coords are:\n";
    std::cout << "longitude: " << 55 << " ~radians~> " << moscow.longitude << "\n"; 
    std::cout << "latitude: " << 37 << " ~radians~> " << moscow.latitude << "\n";
    std::cout << "altitude: " << moscow.altitude << "\n";

    CartesianCoordinates cart_msc = geod_to_cart(moscow);

    std::cout << "Moscow cart coords are:\n";
    std::cout << "x: " << cart_msc.x << "\n"; 
    std::cout << "y: " << cart_msc.y << "\n";
    std::cout << "z: " << cart_msc.z << "\n";

    NewtonsResult res = newton_method(cart_msc);

    std::cout << "Moscow newton res:\n";
    std::cout << "root: " << res.root << "\n"; 
    std::cout << "iters: " << res.iterations << "\n";
    std::cout << "converged: " << res.converged << "\n";

    GeodeticCoordinates geod_msc = cart_to_geod(cart_msc);
    
    std::cout << "Moscow transformed geod coords are:\n";
    std::cout << "longitude: " << geod_msc.longitude*180/M_PI << " ~radians~> " << 
    geod_msc.longitude << (fabs(geod_msc.longitude - moscow.longitude) < 1e-6 ? " (as expected)" : "-") << "\n"; 
    std::cout << "latitude: " << geod_msc.latitude*180/M_PI << " ~radians~> " << 
    geod_msc.latitude << (fabs(geod_msc.latitude - moscow.latitude) < 1e-6 ? " (as expected)" : "-") << "\n"; 
    std::cout << "altitude: " << geod_msc.altitude << "\n";
}