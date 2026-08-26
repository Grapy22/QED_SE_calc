#ifndef UTILS_H_
#define UTILS_H_

#include <iostream>
#include <iomanip>
#include <vector>
#include <quadmath.h>

using std::cout, std::endl, std::cin, std::setprecision, std::pair, std::vector;
using std::make_pair, std::abs;
// --- КОНСТАНТЫ ПРОЕКТА ---

constexpr __float128 alpha = 1.0Q / 137.035999084Q;
constexpr __float128 eps = 1e-14Q;

namespace constants {
    inline constexpr __float128 pi = M_PIq;

    // Без деления
    inline constexpr __float128 c_2pi   = 2 * pi;

    // Делитель 3
    inline constexpr __float128 c_pi_3  = pi / 3.0Q;
    inline constexpr __float128 c_2pi_3 = 2.0Q * pi / 3.0Q;
    inline constexpr __float128 c_4pi_3 = 4.0Q * pi / 3.0Q;

    // Делитель 5
    inline constexpr __float128 c_pi_5  = pi / 5.0Q;
    inline constexpr __float128 c_2pi_5 = 2.0Q * pi / 5.0Q;
    inline constexpr __float128 c_3pi_5 = 3.0Q * pi / 5.0Q;

    // Делитель 7
    inline constexpr __float128 c_pi_7  = pi / 7.0Q;

    // Делитель 15
    inline constexpr __float128 c_pi_15  = pi / 15.0Q;
    inline constexpr __float128 c_2pi_15 = 2.0Q * pi / 15.0Q;
    inline constexpr __float128 c_4pi_15 = 4.0Q * pi / 15.0Q;

    // Делитель 35
    inline constexpr __float128 c_pi_35  = pi / 35.0Q;
    inline constexpr __float128 c_6pi_35 = 6.0Q * pi / 35.0Q;
}


// --- УГЛОВЫЕ ИНТЕГРАЛЫ ---
__float128 Y25(int k_a);
__float128 Y26(int k_a);
__float128 Y34(int k_a);
__float128 Y35(int k_a);
void printq(const __float128& r);

#endif /* UTILS_H_ */
