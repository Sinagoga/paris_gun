#include <iostream>
#include "utils.hpp"
#include "coords.hpp"
#define _USE_MATH_DEFINES
#include <cmath>

using namespace ballistics;
int main() {
    GeodeticCoordinates moscow;
    moscow.latitude = 55*M_PI / 180.0;
    moscow.longitude = 37*M_PI / 180.0;
    moscow.altitude = 124.0;

    std::cout << "Moscow geod coords are:\n";
    std::cout << "longitude: " << moscow.longitude << "\n"; 
    std::cout << "latitude: " << moscow.latitude << "\n";
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
    std::cout << "longitude: " << geod_msc.longitude*180/M_PI << "\n"; 
    std::cout << "latitude: " << geod_msc.latitude*180/M_PI << "\n";
    std::cout << "altitude: " << geod_msc.altitude << "\n";
}