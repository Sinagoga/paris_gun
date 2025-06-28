#ifndef SOLVER_HPP_
#define SOLVER_HPP_

#include <cmath>

namespace ballistics {

/// Параметры для решения ODE
struct BallisticParams {
    double mass;   ///< масса снаряда, кг
    double area;   ///< эффективная площадь поперечного сечения, м^2
};

/// Функция расчёта аэродинамического коэффициента cD(M) = a_drag / (M^2 - 1)^b_drag
double drag_coefficient(double M);

} // namespace ballistics

#endif // SOLVER_HPP_
