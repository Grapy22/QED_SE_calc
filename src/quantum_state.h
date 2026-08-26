#ifndef QUANTUM_STATE_H
#define QUANTUM_STATE_H

#include "utils.h"

using pFunc = __float128 (*)(const __float128&, const __float128&, const __float128&, const __float128&, const __float128&);

class quantum_state
{
private:
    int k, nr, Z, l, l_ ;
    __float128 alphaZ, j , a , N , gamma;
    pair <vector<__float128>, vector<__float128>> C1C2;
protected:
    __float128 en;
public:
    quantum_state(int k , int nr, int Z);
    void ka2jl(const int, int& , int&) ;
    __float128 get_en () const { return en;}
    __float128 get_alphaZ() const {return alphaZ;}
    pair<vector<__float128>, vector<__float128>> C1C2_calc() const;
    vector<__float128> Dm_calc (const __float128& p, const int lmax, const __float128& b = 0.0Q ) const;
    void dirac_pnt_p_dual (const __float128& p, __float128& g, __float128& f,
    const __float128& b = __float128(0), __float128* dg = nullptr, __float128* df = nullptr) const;
    void coeff (const __float128& p, __float128* A = nullptr, __float128* b1 = nullptr, __float128* b2 =nullptr,
    __float128* b3 = nullptr, __float128* a1 = nullptr, __float128* a2 = nullptr, __float128* a3 = nullptr) const;
    const __float128&  ver1 (const __float128&  p) const;
    void prep (const __float128&  x, const __float128&  y, const __float128&  p1, const __float128&  p2,
         const __float128&  ksi, __float128* N = nullptr, __float128* a_ = nullptr ) const;
    void prep_to_1pot_ver1(const __float128&  x, const __float128&  y, __float128&  A0,__float128& C1,
        __float128& C2,__float128& D1, __float128& D2,__float128& F1,__float128& F2,
          __float128&  G1,  __float128&  G2,__float128& H1) const;
    void P_i (const __float128&  x, const __float128&  y, const __float128&  p1, const __float128&  p2,
        const __float128&  ksi,__float128&  g1, __float128&  f1, __float128&  g2, __float128& f2,__float128& P1,
         __float128& P2, __float128& P3, __float128& P4, __float128&  P5, __float128& P6) const;
    void  prep_to_1pot_ver2(const __float128&  x, const __float128&  y, __float128& B, __float128& C,__float128& D,__float128& H) const;
    void R_i (const __float128&  x, const __float128&  y, const __float128&  p1, const __float128&  p2,
        const __float128&  ksi, __float128& g1, __float128& f1, __float128& g2, __float128& f2,__float128&  R1,
             __float128& R2,__float128& R3,__float128& R4, __float128& R5, __float128& R6) const;
    __float128 pot1_ver1 (const __float128&  x, const __float128&  y, const __float128&  p1, const __float128&  p2,
        const __float128&  ksi) const;
    __float128 pot1_ver1_I5 (const __float128&  x, const __float128&  y, const __float128&  p1, const __float128&  p2,
        const __float128&  ksi) const;
    __float128  pot1_ver1_prq (const __float128&  x, const __float128&  y, const __float128&  p, const __float128&  r,
        const __float128&  q) const;
    inline __float128 G (const __float128&  p1, const __float128&  p2) const;
    const __float128& F (const __float128&  x, const __float128&  y, const __float128&  p1, const __float128&  p2,
         const __float128&  ksi,__float128& g1,__float128& f1, __float128&  g2, __float128& f2) const;
    __float128 Y62 (const int, const __float128&  ksi) const;
    __float128 Y63 (const int, const __float128&  ksi) const;
    __float128 Y64 (const int, const __float128&  ksi) const;
    __float128 Y76(const int, const __float128&  ksi, const __float128&  p1, const __float128&  p2) const;
    __float128 Y77(const int, const __float128&  ksi) const;
    __float128 R (const __float128&  x, const __float128&  y, const __float128&  p1, const __float128&  p2,
         const __float128&  ksi,__float128& g1,__float128& f1,__float128& g2,__float128& f2) const;
    __float128 pot1_ver2 (const __float128&  x, const __float128&  y, const __float128&  p1, const __float128&  p2,
        const __float128&  ksi) const;
    __float128 pot1_ver2_prq(const __float128& x, const __float128& y, const __float128& p, const __float128& r,
        const __float128& q) const ;
    __float128 Func_prq(pFunc, const __float128& x, const __float128& y, const __float128& p, const __float128& r,
        const __float128& q) const ;
};



// --- Интегранды для Cuba ---
struct cuba_init5
{
    quantum_state* st;
    __float128  (quantum_state::*func) (const __float128& , const __float128& , const __float128& , const __float128& , const __float128& ) const;
};
struct cuba_init4
{
    quantum_state* st;
    __float128 (quantum_state::*func) (const __float128& , const __float128& , const __float128& , const __float128& ) const;
};
struct cuba_init3
{
    quantum_state* st;
    __float128  (quantum_state::*func) (const __float128& , const __float128& , const __float128& ) const;
};
void Cuba_Int5(cuba_init5 &st_func, const __float128&  epsrel, const __float128&  epsabs ,__float128& res ,__float128& err);
void Cuba_Int3(cuba_init3 &st_func, const __float128&  epsrel, const __float128&  epsabs,__float128& res,__float128& err);
void Cuba_Int4(cuba_init4 &st_func, const __float128&  epsrel, const __float128&  epsabs,__float128& res,__float128& err);
void pot1_ver1(quantum_state& state, const __float128&  abseps, const __float128&  releps,__float128& res,__float128& err);
// void pot1_ver2(quantum_state& state, const __float128&  abseps, const __float128&  releps,__float128& res,__float128& err);
void pot1_ver1_calc(quantum_state& state, const __float128&  abseps,const __float128&  releps);
void pot1_ver2_calc(quantum_state& state,const __float128& abseps,const __float128& releps);
void pot1_calc(quantum_state& state,const __float128& abseps,const __float128& releps);

#endif
