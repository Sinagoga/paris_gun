#include <cmath>
#define _USE_MATH_DEFINES
#include "coords.hpp"
#include "utils.hpp"

namespace ballistics {

CartesianCoordinates geod_to_cart(const GeodeticCoordinates geod_coords) {
    double cos_lat = cos(geod_coords.latitude);
    double sin_lat = sin(geod_coords.latitude);
    double cos_lon = cos(geod_coords.longitude);
    double sin_lon = sin(geod_coords.longitude);
    double n_lat = R_eq/sqrt(1-E_2*sin_lat*sin_lat);

    CartesianCoordinates cart_coords;

    cart_coords.x = (n_lat + geod_coords.altitude)*cos_lat*cos_lon;
    cart_coords.y = (n_lat + geod_coords.altitude)*cos_lat*sin_lon;
    cart_coords.z = (((R_pol*R_pol)/(R_eq*R_eq))*n_lat + geod_coords.altitude)*sin_lat;

    return cart_coords;
};

GeodeticCoordinates cart_to_geod(const CartesianCoordinates cart_coords) {
    NewtonsResult newton_res = newton_method(cart_coords);
    double t = newton_res.root;
    double lat = atan(t);

    GeodeticCoordinates geod_coords;

    geod_coords.longitude = atan2(cart_coords.y, cart_coords.x);
    geod_coords.latitude = lat;
    geod_coords.altitude = calculate_h(cart_coords, lat);

    return geod_coords;
};

}
