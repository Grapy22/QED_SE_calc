
#include "quantum_state.h"
extern "C" {
    // 1. Создаем правильный тип для указателя на 128-битную функцию
    typedef int (*quad_integrand_t)(const int *ndim, const __float128 *x, const int *ncomp, __float128 *f, void *userdata);
    
    // 2. Объявляем саму функцию Cuhre, где все типы данных заменены на __float128
    void Cuhre(const int ndim, const int ncomp,
               quad_integrand_t integrand, void *userdata,
               const int nvec, const __float128 epsrel, const __float128 epsabs,
               const int flags, const int mineval, const int maxeval,
               const int key, const char *statefile, void *spin,
               int *nregions, int *neval, int *fail,
               __float128 integral[], __float128 error[], __float128 prob[]);
}
int Cuba_integrand5_to_p0(const int *ndim, const __float128 x[],
                          const int *ncomp, __float128 F[], void *userdata)
{
    auto *st_func = static_cast<cuba_init5*>(userdata);
    quantum_state* st = st_func->st;
    __float128 p0 = 2 * st->get_alphaZ();
    __float128 p = p0 * x[0]; // p in [0,p0]
    __float128 r = p * x[1];
    __float128 abs_r = fabsq(r);                               // r in [-p, p]
    __float128 q = sqrtq(2.0Q) * (abs_r + (p - abs_r) * x[2]); // q in [√2|r|, √2p]

    __float128 x_var = x[3];
    __float128 y_var = x[4];
    __float128 jac = p0 * p * (sqrtq(2.0Q) * (p - fabsq(r)));
    F[0] = (st->*(st_func->func))(x_var, y_var, p, r, q) * jac;

    r = -r;

    F[0] += (st->*(st_func->func))(x_var, y_var, p, r, q) * jac;

    return 0;
}

int Cuba_integrand5_to_inf(const int *ndim, const __float128 x[],
                           const int *ncomp, __float128 F[], void *userdata)
{
    auto *st_func = static_cast<cuba_init5*>(userdata);
    quantum_state* st = st_func->st;
    __float128 p0 = 2 * st->get_alphaZ();
    __float128 p = p0 / x[0];

    __float128 r = p * x[1]; // r in [-p, p]
    __float128 abs_r = fabsq(r);
    __float128 q = sqrtq(2.0Q) * (abs_r + (p - abs_r) * x[2]); // q in [√2|r|, √2p]
    __float128 x_var = x[3];
    __float128 y_var = x[4];
    __float128 jac = (p * p / p0) * p * sqrtq(2.0Q) * (p - abs_r);
    F[0] = (st->*(st_func->func))(x_var, y_var, p, r, q) * jac;

    r = -r;

    F[0] += (st->*(st_func->func))(x_var, y_var, p, r, q) * jac;

    return 0;
}

void Cuba_Int5(cuba_init5 &st_func, const __float128& epsrel,const __float128& epsabs, __float128 &res, __float128 &err)
{
    const int ndim = 5;  // Размерность
    const int ncomp = 1; // Количество компонент результата

    const int flags = 0;           // Флаги (0 - стандартные настройки)
    const int mineval = 1;         // Минимальное число вычислений
    const int maxeval = 256000000; // Максимальное число вычислений
    const int key = 0;             // Алгоритм интегрирования

    // Выходные переменные
    int nregions1, neval1, fail1;
    __float128 integral1[1], error1[1], prob1[1];
    int nregions2, neval2, fail2;
    __float128 integral2[1], error2[1], prob2[1];
    Cuhre(ndim, ncomp, Cuba_integrand5_to_p0, &st_func, 1, // 1 - nvec (число векторов)
          epsrel, epsabs, flags,
          mineval, maxeval, key,
          NULL, NULL, // statefile, spin (не используем)
          &nregions1, &neval1, &fail1,
          integral1, error1, prob1);
    __float128 res1 = integral1[0];
    // __float128  fctr = (4 * alpha) / pow(2 * pi, 6) / pi;
    // cout << setprecision(15)<< "[0 ; p0] :" << res1 * fctr<< '\n';
    __float128 err1 = error1[0];
    Cuhre(ndim, ncomp, Cuba_integrand5_to_inf, &st_func, 1, // 1 - nvec (число векторов)
          epsrel, epsabs, flags,
          mineval, maxeval, key,
          NULL, NULL, // statefile, spin (не используем)
          &nregions2, &neval2, &fail2,
          integral2, error2, prob2);
    __float128 res2 = integral2[0];
    __float128 err2 = error2[0];
    // cout << setprecision(15)<< "[p0 ; inf] :" << res2 * fctr << '\n';

    res = res1 + res2;
    err = sqrtq(err1 * err1 + err2 * err2);
}
                          
// int Cuba_integrand3_to_p0(const int *ndim, const __float128 x[],
//                           const int *ncomp, __float128 F[], void *userdata)
// {
//     auto *st_func = static_cast<cuba_init3*>(userdata);
//     quantum_state* st = st_func->st;
//     __float128 p0 = 2 * st->get_alphaZ();
//     __float128 x_var = x[0];
//     __float128 y_var = x[1];
//     __float128 p = p0 * x[2];

//     F[0] = (st->*(st_func->func))(x_var, y_var, p) * p0;
//     return 0;
// }
// int Cuba_integrand3_to_inf(const int *ndim, const __float128 x[],
//                           const int *ncomp, __float128 F[], void *userdata)
// {
//     auto *st_func = static_cast<cuba_init3*>(userdata);
//     quantum_state* st = st_func->st;
//     __float128 p0 = 2 * st->get_alphaZ();
//     __float128 x_var = x[0];
//     __float128 y_var = x[1];
//     __float128 p = p0 / x[2];

//     F[0] = (st->*(st_func->func))(x_var, y_var, p) * p * p / p0;
//     return 0;
// }

// void Cuba_Int3(cuba_init3 &st_func,const __float128& epsrel,const __float128& epsabs, __float128 &res, __float128 &err)
// {
//     const int ndim = 3;  // Размерность
//     const int ncomp = 1; // Количество компонент результата

//     const int flags = 0;           // Флаги (0 - стандартные настройки)
//     const int mineval = 1;         // Минимальное число вычислений
//     const int maxeval = 256000000; // Максимальное число вычислений
//     const int key = 0;             // Алгоритм интегрирования

//     // Выходные переменные
//     int nregions1, neval1, fail1;
//     __float128 integral1[1], error1[1], prob1[1];
//     int nregions2, neval2, fail2;
//     __float128 integral2[1], error2[1], prob2[1];
//     Cuhre(ndim, ncomp, Cuba_integrand3_to_p0, &st_func, 1, 
//           epsrel, epsabs, flags,
//           mineval, maxeval, key,
//           NULL, NULL, 
//           &nregions1, &neval1, &fail1,
//           integral1, error1, prob1);
//     __float128 res1 = integral1[0];
//     __float128 err1 = error1[0];   
//     Cuhre(ndim, ncomp, Cuba_integrand3_to_inf, &st_func, 1, 
//           epsrel, epsabs, flags,
//           mineval, maxeval, key,
//           NULL, NULL, // statefile, spin (не используем)
//           &nregions2, &neval2, &fail2,
//           integral2, error2, prob2);
//     __float128 res2 = integral2[0];
//     __float128 err2 = error2[0];
//     res = res1 + res2;
//     err = sqrt(err1 * err1 + err2 * err2);
// }
// int Cuba_integrand4_to_p0(const int *ndim, const __float128 x[],
//                           const int *ncomp, __float128 F[], void *userdata)
// {
//     auto *st_func = static_cast<cuba_init4*>(userdata);
//     quantum_state* st = st_func->st;
//     __float128 p0 = 2 * st->get_alphaZ();
//     __float128 x_var = x[0];
//     __float128 y_var = x[1];
//     __float128 u = p0 * x[2];
//     __float128 t = x[3];

//     F[0] = (st->*(st_func->func))(x_var, y_var,  u, t) * p0; // jacobian = p0

//     t = -t;

//     F[0] += (st->*(st_func->func))(x_var, y_var, u, t) * p0;
//     return 0;
// }
// int Cuba_integrand4_to_inf(const int *ndim, const __float128 x[],
//                           const int *ncomp, __float128 F[], void *userdata)
// {
//     auto *st_func = static_cast<cuba_init4*>(userdata);
//     quantum_state* st = st_func->st;
//     __float128 p0 = 2 * st->get_alphaZ();
//     __float128 x_var = x[0];
//     __float128 y_var = x[1];
//     __float128 u = p0 / x[2];
//     __float128 t = x[3];

//     F[0] = (st->*(st_func->func))(x_var, y_var,  u, t) * u * u / p0;

//     t = -t;

//     F[0] += (st->*(st_func->func))(x_var, y_var,  u, t) * u * u / p0;
//     return 0;
// }

// void Cuba_Int4(cuba_init4 &st_func,con __float128 epsrel, __float128 epsabs, __float128 &res, __float128 &err)
// {
//     const int ndim = 4;  // Размерность
//     const int ncomp = 1; // Количество компонент результата

//     const int flags = 0;           // Флаги (0 - стандартные настройки)
//     const int mineval = 1;         // Минимальное число вычислений
//     const int maxeval = 256000000; // Максимальное число вычислений
//     const int key = 0;             // Алгоритм интегрирования

//     // Выходные переменные
//     int nregions1, neval1, fail1;
//     __float128 integral1[1], error1[1], prob1[1];
//     int nregions2, neval2, fail2;
//     __float128 integral2[1], error2[1], prob2[1];
//     Cuhre(ndim, ncomp, Cuba_integrand4_to_p0, &st_func, 1, 
//           epsrel, epsabs, flags,
//           mineval, maxeval, key,
//           NULL, NULL, 
//           &nregions1, &neval1, &fail1,
//           integral1, error1, prob1);
//     __float128 res1 = integral1[0];
//     __float128 err1 = error1[0];   
//     Cuhreq(ndim, ncomp, Cuba_integrand4_to_inf, &st_func, 1, 
//           epsrel, epsabs, flags,
//           mineval, maxeval, key,
//           NULL, NULL, // statefile, spin (не используем)
//           &nregions2, &neval2, &fail2,
//           integral2, error2, prob2);
//     __float128 res2 = integral2[0];
//     __float128 err2 = error2[0];
//     res = res1 + res2;
//     err = sqrtq(err1 * err1 + err2 * err2);
// }