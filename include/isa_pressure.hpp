#ifndef ISA_PRESSURE_HPP_
#define ISA_PRESSURE_HPP_

#include <cmath>

namespace ballistics {

/**
 * @brief International Standard Atmosphere (ISA) pressure model
 *
 *  - 0 m ≤ h ≤ 11 000 m:
 *      p(h) = p0 * (T(h)/T0)^(g0/(R*L))
 *  - 11 000 m < h ≤ 20 000 m:
 *      p(h) = p_trop * exp(-g0 * (h - h_trop) / (R * T_trop))
 */
class ISAPressure {
public:
    // Note: not constexpr due to use of std::pow
    static double at(double h) {
        const double p0     = 101325.0;    // Pa at sea level
        const double T0     = 288.15;      // K at sea level
        const double L      = 0.0065;      // K/m lapse rate
        const double h_trop = 11000.0;     // m tropopause height
        const double R      = 287.05;      // J/(kg·K)
        const double g0     = 9.80665;     // m/s² gravity
        // Exponent for barometric formula in troposphere
        const double expn   = g0 / (R * L);

        const double T_trop = T0 - L * h_trop;
        const double p_trop = p0 * std::pow(T_trop / T0, expn);

        if (h < 0.0) {
            return p0;
        } else if (h <= h_trop) {
            double T = T0 - L * h;
            return p0 * std::pow(T / T0, expn);
        } else {
            return p_trop * std::exp(-g0 * (h - h_trop) / (R * T_trop));
        }
    }
};

} // namespace ballistics

#endif // ISA_PRESSURE_HPP_
