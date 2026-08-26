#include "quantum_state.h"


void quantum_state::prep (const __float128&  x, const __float128& y,const __float128&  p1,
    const __float128&  p2, const __float128& ksi, __float128 * N , __float128 * a_  ) const
{
    __float128  pp_1 = en * en - p1 * p1 ;
    __float128  pp_2 = en * en - p2 * p2 ;
    __float128  ro1 = 1.0Q - pp_1;
    __float128  ro2 = 1.0Q - pp_2;
    __float128  b2 = y * y * pp_1 + (1.0Q - y) * (1.0Q - y) * pp_2 + 2 * y * (1.0Q - y) * (en * en - p1 * p2 * ksi);

    if (N){
        *N = x * b2 + y * ro1 + (1.0Q - y) * ro2;
    }
    if (a_){
        *a_ = x * y * (1.0Q - x * y) * pp_1 + x * (1.0Q - y) * (1.0Q - x + x * y) * pp_2 
        - 2 * (1.0Q - x * y) * (1.0Q - x + x * y) * (en * en - p1 * p2 * ksi);
    }
}         

void quantum_state::prep_to_1pot_ver1 (const __float128& x, const __float128& y, __float128& A0, __float128& C1,
     __float128& C2, __float128& D1, __float128& D2, __float128& F1, __float128& F2, __float128& G1, __float128& G2, __float128& H1) const
{
    __float128 m_1x = 1.0Q - x ;
    __float128 m_1xy = 1.0Q - x * y;
    __float128 m_bracket = 1.0Q - x + x * y;
    __float128 m_x1y = x * (1.0Q - y);  

    A0 =  1.0Q + 2 * en * en * m_1x * m_1xy;
    C1 = x * y * m_1xy;
    C2 = -x * x * y * (1.0Q - y); 
    D1 = -m_bracket * m_1xy;
    D2 = m_bracket * m_x1y;
    F1 = -en * m_1x * m_1xy;
    F2 = -en * m_1x * m_x1y;
    G1 = 2.0Q;
    G2 = -G1;
    H1 = -2 * en * m_1x * m_1xy;
}

void quantum_state::P_i (const __float128& x, const __float128& y, const __float128& p1, const __float128& p2,
    const __float128& ksi,__float128& g1, __float128& f1, __float128& g2, __float128& f2, __float128& P1,
        __float128& P2, __float128& P3, __float128& P4, __float128& P5, __float128& P6) const
{
    __float128 N, A0, C1, C2, D1, D2, F1, F2, G1, G2, H1;
    prep(x , y, p1, p2, ksi, &N);
    prep_to_1pot_ver1(x, y, A0, C1, C2, D1, D2, F1, F2, G1, G2, H1);

    __float128 fctr = (1.0Q - y) / (N * N);

    P1 = fctr * (A0 * g1 * g2 + F1 * (en * g1 + p1 * f1)*g2 + F2 * g1 * (en * g2 + p2 * f2));

    P2 = fctr * (A0 * f1 * f2 + F1 * (en * f1 + p1 * g1)*f2 + F2 * f1 * (en * f2 + p2 * g2));

    P3 = fctr * (C1 * (en * g1 + p1 * f1)*f2 + C2 * g1 * (en * f2 + p2 * g2) + (H1 - G1) * g1 *f2);

    P4 = fctr * (C1 * (en * f1 + p1 * g1)*g2 + C2 * f1 * (en * g2 + p2 * f2) + (H1 + G1) * f1 * g2);

    P5 = fctr * (D1 * (en * g1 + p1 * f1)*f2 + D2 * g1 * (en * f2 + p2 * g2) - G2 * g1 * f2);

    P6 = fctr * (D1 * (en * f1 + p1 * g1)*g2 + D2 * f1 * (en * g2 + p2 * f2) + G2 * f1 * g2);
}        
void  quantum_state::prep_to_1pot_ver2(const __float128& x, const __float128& y,__float128& B,  __float128& C,  __float128& D,  __float128& H) const 
{
    __float128 m_1x = 1.0Q - x ;
    __float128 m_1xy = 1.0Q - x * y;
    __float128 m_bracket = 1.0Q - x + x * y;
    B = 2.0Q * m_1xy * m_1x ;
    C = 2.0Q * m_bracket * m_1x ;
    D = -m_1x ;
    H = -4.0Q * m_1x ;
}
void quantum_state::R_i (const __float128& x, const __float128& y, const __float128& p1, const __float128& p2,
    const __float128& ksi,  __float128& g1, __float128& f1, __float128& g2, __float128& f2, __float128& R1,
        __float128& R2, __float128& R3, __float128& R4, __float128& R5, __float128& R6) const
{
    __float128 a_, N, B, C, D, H;
    prep(x , y, p1, p2, ksi, &N, &a_);
    prep_to_1pot_ver2(x, y, B, C, D, H);
    
    __float128 K1, K2;
    K1 = (a_ + 1.0Q + 2 * x * N + en * H)*g1*g2 + en * B * (en * g1 + p1 * f1)*g2 
    + en * C * g1 * (en * g2 + p2 * f2) + D * (en * g1 + p1 * f1) * (en * g2 + p2 * f2) ;
    K2 = (a_ + 1.0Q + 2 * x * N - en * H)*f1*f2 + en * B * (en * f1 + p1 * g1)*f2 
    + en * C * f1 * (en * f2 + p2 * g2) + D * (en * f1 +p1 * g1) * (en * f2 + p2 * g2);

    __float128 fctr = 1.0Q / N ;
    __float128 fctr2 = 2 * (1.0Q - y) * fctr * fctr;

    R1 = fctr * (en * C * g1 * f2 + D * (en * g1 + p1 * f1)*f2);

    R2 = fctr * (en * C * f1 * g2 + D * (en * f1 + p1 * g1)*g2);

    R3 =  -fctr2 * x * y * K1 - 2 * (1.0Q - x) * g1 * g2 * fctr;

    R4 = fctr2 * (1.0Q - x + x * y) * K1;

    R5 = -fctr2 * x * y * K2 - 2 * (1.0Q - x) * fctr * f1 * f2;

    R6 = fctr2 * (1.0Q - x + x * y) * K2;
}        
