
#include "quantum_state.h"

__float128 quantum_state::pot1_ver1(const __float128& x, const __float128& y, const __float128& p1, const __float128& p2,
        const __float128& ksi) const
{
    __float128 g1, f1, g2, f2;
    dirac_pnt_p_dual(p1, g1, f1);
    dirac_pnt_p_dual(p2, g2, f2);

    __float128 P1 , P2, P3, P4 , P5, P6;
    P_i(x , y , p1 , p2 , ksi , g1 , f1 , g2, f2, P1 , P2, P3, P4 , P5, P6);

    __float128 q2 = p1 * p1 + p2 * p2 - 2 * p1 * p2 * ksi;

    __float128 result;
    result = P1 * Y62(k , ksi) + P2 * Y62(-k , ksi) + p1 * ( P3 * Y63(k , ksi)
    + P4 * Y63(-k , ksi)) + p2 * (P5 * Y64(k , ksi) + P6 * Y64(-k , ksi));
    return  result / q2;
}

__float128 quantum_state::pot1_ver1_prq(const __float128& x, const __float128& y, const __float128& p, const __float128& r,
        const __float128& q) const // сделать аргумент как указатель на функцию
{
    __float128 p1 = (p + r) / M_SQRT2q;
    __float128 p2 = (p - r) / M_SQRT2q;
    __float128 ksi = (p1*p1 + p2*p2 - q*q) / (2 * p1 * p2);
    return pot1_ver1(x, y, p1, p2, ksi) * q ;
}



// inline double quantum_state::G (double p1, double p2) const
// {
//     return log((p1 + p2) / std::abs(p1 - p2)) ;
// }

// double quantum_state::F (double x, double y, double p1, double p2, double ksi, double& g1,
//           double& f1,  double& g2, double& f2) const
// {
//     double P1 , P2, P3, P4 , P5, P6;
//     P_i(x , y , p1 , p2 , ksi , g1 , f1 , g2, f2, P1 , P2, P3, P4 , P5, P6);
//     double res;
//     res = P1 * Y62(k , ksi) + P2 * Y62(-k , ksi) + p1 * ( P3 * Y63(k , ksi)
//     + P4 * Y63(-k , ksi)) + p2 * (P5 * Y64(k , ksi) + P6 * Y64(-k , ksi));
//     return res;
// }
// double quantum_state::pot1_ver1_I5(double x, double y, double p1, double p2,
//         double ksi) const
// {
//     double g1, f1, g2, f2;
//     dirac_pnt_p_dual(p1, g1, f1);
//     dirac_pnt_p_dual(p2, g2, f2);

//     double q2 = p1 * p1 + p2 * p2 - 2 * p1 * p2 * ksi;
//     double Vq = -4 * pi * alphaZ / q2;

//     double result = F (x, y, p1, p2, ksi, g1, f1, g2, f2);
//     double result_s = F (x, y, p1, p2, 1.0, g1, f1, g2, f2);

//     return  (result - result_s) * Vq ;
// }
// double quantum_state::pot1_ver1_I4(double x, double y, double p1,
//     double p2) const
// {
//     double g1, f1, g2, f2;
//     dirac_pnt_p_dual(p1, g1, f1);
//     dirac_pnt_p_dual(p2, g2, f2);

//     double result = F(x, y, p1, p2, 1.0, g1, f1, g2, f2) ;
//    // double result_s = F(x, y, p1, p1, 1.0, g1, f1, g1, f1) ;
//     return  result * G(p1 , p2) ;
// }

// double quantum_state::pot1_ver1_I3(double x, double y, double p1) const
// {
//     double g1, f1;
//     dirac_pnt_p_dual(p1, g1, f1);
//     double result = F(x, y, p1, p1, 1.0, g1, f1, g1, f1) ;
//    //double G_an = pi * pi / (2 * p1) ;
//     return  result * p1;
// }
