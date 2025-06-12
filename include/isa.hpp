#ifndef ISA_HPP_
#define ISA_HPP_

#include <cmath>

namespace ballistics {

/**
 * @brief International Standard Atmosphere (ISA) temperature model
 *
 *  - 0 m ≤ h ≤ 11 000 m: T(h) = T0 − L·h
 *  - 11 000 m < h ≤ 20 000 m: T(h) = T_tropopause (постоянная)
 *
 * Всё — constexpr, можно использовать в compile‑time.
 */

class ISA {
    public:
    static constexpr double get_temperature(double h) {
        constexpr double T0     = 288.15;    // K, у уровня моря
        constexpr double L      = 0.0065;    // K/m, градиент в тропосфере
        constexpr double h_trop = 11000.0;   // m, высота тропопаузы
        constexpr double T_trop = T0 - L * h_trop;

        return (h < 0.0)         ? T0
                : (h <= h_trop)     ? (T0 - L * h)
                : /* h ≤ 20 000 */   T_trop;
    };

    static constexpr double get_pressure(double h) {
        const double p0     = 101325.0;    // Pa at sea level
        const double T0     = 288.15;      // K at sea level
        const double L      = 0.0065;      // K/m lapse rate
        const double h_trop = 11000.0;     // m tropopause height
        const double R      = 287.0528;    // J/(kg·K)
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
    };

    static constexpr double get_density(double h) {
        constexpr double R_air = 287.0528; // J/(kg·K) specific gas constant for the air

        return ISA::get_pressure(h) / (R_air * ISA::get_temperature(h));
    };
};

};

#endif