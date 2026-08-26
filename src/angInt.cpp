#include "quantum_state.h"

__float128 Y25(int k_a) {
    // Сдвигаем индекс, чтобы k_a от -3 до 3 превратился в индекс от 0 до 6
    if (k_a < -3 || k_a > 3 || k_a == 0) {
        std::cerr << "Error: k_a out of range in Y25 (k_a = " << k_a << ")" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Значения для k_a = [-3, -2, -1,  0,  1,  2,  3]
    static constexpr __float128 lookup[] = {
         1.0Q / 5.0Q,  // -3
         1.0Q / 3.0Q,  // -2
         1.0Q,         // -1
         0.0Q,         //  0 (не используется)
        -1.0Q / 3.0Q,  //  1
        -1.0Q / 5.0Q,  //  2
        -1.0Q / 7.0Q   //  3
    };

    return lookup[k_a + 3];
}

__float128 Y26(int k_a) {
    if (k_a < -3 || k_a > 3 || k_a == 0) {
        std::cerr << "Error: k_a out of range in Y26 (k_a = " << k_a << ")" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    // Индексы: [-3, -2, -1, 0, 1, 2, 3] -> сдвиг +3
    static constexpr __float128 lookup[] = {
        -6.0Q / 35.0Q,  // -3
        -4.0Q / 15.0Q,  // -2
        -2.0Q / 3.0Q,   // -1
         0.0Q,          //  0 (не используется)
         2.0Q / 3.0Q,   //  1
         4.0Q / 15.0Q,  //  2
         6.0Q / 35.0Q   //  3
    };
    return lookup[k_a + 3];
}

__float128 Y34(int k_a) {
    if (k_a < -3 || k_a > 3 || k_a == 0) {
        std::cerr << "Error: k_a out of range in Y34 (k_a = " << k_a << ")" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    static constexpr __float128 lookup[] = {
        -2.0Q / 5.0Q,   // -3
         1.0Q / 3.0Q,   // -2
         0.0Q,          // -1
         0.0Q,          //  0 (не используется)
         2.0Q / 3.0Q,   //  1
         3.0Q / 5.0Q,   //  2
         4.0Q / 7.0Q    //  3
    };
    return lookup[k_a + 3];
}

__float128 Y35(int k_a) {
    if (k_a < -3 || k_a > 3 || k_a == 0) {
        std::cerr << "Error: k_a out of range in Y35 (k_a = " << k_a << ")" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    static constexpr __float128 lookup[] = {
        -24.0Q / 35.0Q, // -3
        -4.0Q / 5.0Q,   // -2
        -4.0Q / 3.0Q,   // -1
         0.0Q,          //  0 (не используется)
         0.0Q,          //  1
        -4.0Q / 15.0Q,  //  2
        -12.0Q / 35.0Q  //  3
    };
    return lookup[k_a + 3];
}

__float128 quantum_state::Y62(const int k_a,const __float128& ksi) const {
    switch (k_a) {
        case -1: return constants::c_2pi;
        case 1:  return -constants::c_2pi_3 * ksi;
        case -2: return constants::c_2pi_3 * ksi;
        case 2:  return constants::c_pi_5 * (1.0Q - 3.0Q * ksi * ksi);
        case -3: return -constants::c_pi_5 * (1.0Q - 3.0Q * ksi * ksi);
        case 3:  return constants::c_pi_7 * (3.0Q * ksi - 5.0Q * ksi * ksi * ksi);
        default:
            std::cerr << "Error: k_a out of range in Y62 (k_a = " << k_a << ")" << std::endl;
            std::exit(EXIT_FAILURE);
    }
}

__float128 quantum_state::Y63(const int k_a,const __float128& ksi) const {
    switch (k_a) {
        case -1: return -constants::c_2pi_3 * ksi;
        case 1:  return -constants::c_2pi_3;
        case -2: return constants::c_pi_15 * (1.0Q - 3.0Q * ksi * ksi);
        case 2:  return -constants::c_2pi_15 * ksi;
        case -3: return -constants::c_pi_35 * (-3.0Q * ksi + 5.0Q * ksi * ksi * ksi);
        case 3:  return -constants::c_pi_35 * (-1.0Q + 3.0Q * ksi * ksi);
        default:
            std::cerr << "Error: k_a out of range in Y63 (k_a = " << k_a << ")" << std::endl;
            std::exit(EXIT_FAILURE);
    }
}

__float128 quantum_state::Y64(const int k_a, const __float128& ksi) const {
    switch (k_a) {
        case -1: return -constants::c_2pi_3;
        case 1:  return -constants::c_2pi_3 * ksi;
        case -2: return -constants::c_2pi_15 * ksi;
        case 2:  return constants::c_pi_15 * (1.0Q - 3.0Q * ksi * ksi);
        case -3: return constants::c_pi_35 * (1.0Q - 3.0Q * ksi * ksi);
        case 3:  return -constants::c_pi_35 * (-3.0Q * ksi + 5.0Q * ksi * ksi * ksi);
        default:
            std::cerr << "Error: k_a out of range in Y64 (k_a = " << k_a << ")" << std::endl;
            std::exit(EXIT_FAILURE);
    }
}

__float128 quantum_state::Y76(const int k_a, const __float128& ksi,const __float128& p1, const __float128& p2) const {
    switch (k_a) {
        case -1: return -constants::c_4pi_3 * (p2 - p1 * ksi);
        case 1:  return -constants::c_4pi_3 * (p1 - p2 * ksi);
        case -2: return -constants::c_4pi_15 * (p1 + 2.0Q * p2 * ksi - 3.0Q * p1 * ksi * ksi);
        case 2:  return -constants::c_4pi_15 * (p2 + 2.0Q * p1 * ksi - 3.0Q * p2 * ksi * ksi);
        case -3: return constants::c_6pi_35 * (p2 - 3.0Q * p2 * ksi * ksi + p1 * ksi * (-3.0Q + 5.0Q * ksi * ksi));
        case 3:  return constants::c_6pi_35 * (p1 - 3.0Q * p1 * ksi * ksi + p2 * ksi * (-3.0Q + 5.0Q * ksi * ksi));
        default:
            std::cerr << "Error: k_a out of range in Y76 (k_a = " << k_a << ")" << std::endl;
            std::exit(EXIT_FAILURE);
    }
}

__float128 quantum_state::Y77(const int k_a, const __float128& ksi) const {
    switch (k_a) {
        case -1: return 0.0Q;
        case 1:  return constants::c_2pi_3 * (1.0Q - ksi * ksi);
        case -2: return constants::c_pi_3 * (1.0Q - ksi * ksi);
        case 2:  return constants::c_3pi_5 * (1.0Q - ksi * ksi) * ksi;
        case -3: return constants::c_2pi_5 * (1.0Q - ksi * ksi) * ksi;
        case 3:  return -constants::c_pi_7 * (1.0Q - 6.0Q * ksi * ksi + 5.0Q * ksi * ksi * ksi * ksi);
        default:
            std::cerr << "Error: k_a out of range in Y77 (k_a = " << k_a << ")" << std::endl;
            std::exit(EXIT_FAILURE);
    }
}
