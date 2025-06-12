#ifndef COORDS_HPP_
#define COORDS_HPP_

namespace ballistics {

#define F (double)1/298
#define R_eq (double)6378136.3
#define R_pol R_eq*(1 - F)
#define E_2 (R_eq*R_eq - R_pol*R_pol)/(R_eq*R_eq)

struct GeodeticCoordinates {
    double longitude;
    double latitude;
    double altitude;
};

struct CartesianCoordinates {
    double x;
    double y;
    double z;
};

CartesianCoordinates geod_to_cart(const GeodeticCoordinates geod_coords);
GeodeticCoordinates cart_to_geod(const CartesianCoordinates cart_coords);

}

#endif