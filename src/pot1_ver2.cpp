
#include "quantum_state.h"

__float128 quantum_state::R(const __float128& x, const __float128& y, const __float128& p1, const __float128& p2,
         const __float128& ksi,__float128& g1, __float128& f1, __float128& g2,  __float128& f2) const  
{
    __float128 R1, R2, R3, R4, R5, R6;
    R_i (x, y, p1, p2, ksi, g1, f1, g2, f2, R1, R2, R3, R4, R5, R6) ;
    __float128 res = R1 * Y76(k , ksi , p1 , p2) + R2 * Y76(-k , ksi , p1 , p2) +
    p1 * p2 * ( Y77(k , ksi) * (R3 + R4) +  Y77(-k , ksi) * (R5 + R6) );
    return res;
}
        
__float128 quantum_state::pot1_ver2 (const __float128& x, const __float128& y, const __float128& p1, const __float128& p2, const __float128& ksi) const
{
    __float128 g1, f1, g2, f2;
    dirac_pnt_p_dual(p1, g1, f1);
    dirac_pnt_p_dual(p2, g2, f2);
    __float128 q2 = p1 * p1 + p2 * p2 - 2 * p1 * p2 * ksi;
    __float128 res = R(x, y, p1, p2, ksi, g1, f1, g2, f2);
    return res / (q2 * q2);
}   
__float128 quantum_state::pot1_ver2_prq(const __float128& x, const __float128& y, const __float128& p, const __float128& r,
        const __float128& q) const // сделать аргумент как указатель на функцию
{
    __float128 p1 = (p + r) / M_SQRT2q;
    __float128 p2 = (p - r) / M_SQRT2q;
    __float128 ksi = (p1*p1 + p2*p2 - q*q) / (2 * p1 * p2);
    return pot1_ver2(x, y, p1, p2, ksi) * q ;
}
