#ifndef ISA_TEMPERATURE_HPP_
#define ISA_TEMPERATURE_HPP_

namespace ballistics {

/**
 * @brief International Standard Atmosphere (ISA) temperature model
 *
 *  - 0 m ≤ h ≤ 11 000 m: T(h) = T0 − L·h
 *  - 11 000 m < h ≤ 20 000 m: T(h) = T_tropopause (постоянная)
 *
 * Всё — constexpr, можно использовать в compile‑time.
 */
class ISATemperature {
public:
    static constexpr double at(double h) {
        constexpr double T0     = 288.15;    // K, у уровня моря
        constexpr double L      = 0.0065;    // K/m, градиент в тропосфере
        constexpr double h_trop = 11000.0;   // m, высота тропопаузы
        constexpr double T_trop = T0 - L * h_trop;

        return (h < 0.0)         ? T0
             : (h <= h_trop)     ? (T0 - L * h)
             : /* h ≤ 20 000 */   T_trop;
    }
};

} // namespace ballistics

#endif // ISA_TEMPERATURE_HPP_
