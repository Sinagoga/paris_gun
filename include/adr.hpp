#ifndef ADR_HPP_
#define ADR_HPP_

namespace ballistics {

    double c_d(double m, double a, double b);
    double s(double diameter);
    double R(double dens, double v, double s, double c_d);

}

#endif