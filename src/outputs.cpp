
#include "quantum_state.h"

void printq(const __float128& r)
{
    char buf[128];
    quadmath_snprintf(buf, sizeof(buf), "%-#*.30Qf", 40, r);
    printf("%s\n", buf);
}
void print_result(std::ostream& out,
                  const std::string& label,
                  __float128 res,
                  __float128 err,
                  __float128 epsabs,
                  int scale = 0)
{
    // Масштаб
    __float128 scale_factor = powq(10.0Q, (__float128)scale);

    // Аналог Fortran:
    // actual_err = abs(err * scale)
    // scaled_res = res * scale
    __float128 actual_err = fabsq(err * scale_factor);
    __float128 scaled_res = res * scale_factor;

    // eff_eps = max(actual_err, abs(epsabs * scale))
    __float128 eff_eps = fmaxq(actual_err,
                               fabsq(epsabs * scale_factor));

    int digits;

    if (eff_eps > 1e-20Q)
    {
        digits = std::max(0,
                 (int)(-log10q(eff_eps)) + 1);
    }
    else
    {
        digits = 15;
    }

    digits = std::min(digits, 15);

    // Формат вывода
    std::string fmt = "%." + std::to_string(digits) + "Qf";

    char res_buf[256];
    char err_buf[256];

    quadmath_snprintf(res_buf,
                      sizeof(res_buf),
                      fmt.c_str(),
                      scaled_res);

    quadmath_snprintf(err_buf,
                      sizeof(err_buf),
                      fmt.c_str(),
                      actual_err);

    // Вывод
    out << label << " : "
        << res_buf
        << " ± "
        << err_buf;

    if (scale != 0)
        out << " x10^" << scale;

    out << '\n';

    out.flush();
}
void pot1_ver1_calc(quantum_state& state,const __float128& abseps,const __float128& releps)
{
    cuba_init5 data;
    data.st = &state;
    data.func = &quantum_state::pot1_ver1_prq;
    __float128 res, err, fctr;
    fctr = (4 * alpha) / powq(2 * constants::pi, 6) / constants::pi;
    __float128 abseps_f = abseps / fctr;
    fctr = fctr * (-4 * constants::pi * state.get_alphaZ());
    Cuba_Int5(data, releps, abseps_f, res, err);
    print_result(cout, "po1_ver1", res * fctr, err * fctr, abseps, 6);
}
void pot1_ver2_calc(quantum_state& state,const __float128& abseps,const __float128& releps)
{
    cuba_init5 data;
    data.st = &state;
    data.func = &quantum_state::pot1_ver2_prq;
    __float128 res, err, fctr;
    fctr = (4 * alpha) / powq(2 * constants::pi, 7);
    __float128 abseps_f = abseps / fctr;
    fctr = fctr * (8 * constants::pi * state.get_alphaZ());
    Cuba_Int5(data, releps, abseps_f, res, err);
    print_result(cout, "po1_ver2", res * fctr, err * fctr, abseps, 6);
}
void pot1_calc(quantum_state& state,const __float128& abseps,const __float128& releps)
{
    cuba_init5 data;
    data.st = &state;
    data.func = &quantum_state::pot1_ver1_prq;
    __float128 res1, err1, fctr1, res2, err2, fctr2;
    fctr1 = (4 * alpha) / powq(2 * constants::pi, 6) / constants::pi;
    fctr2 = (4 * alpha) / powq(2 * constants::pi, 7);
    __float128 abseps_f1 = abseps / fctr1;
    fctr1 = fctr1 * (-4 * constants::pi * state.get_alphaZ());
    Cuba_Int5(data, releps, abseps_f1, res1, err1);

    data.func = &quantum_state::pot1_ver2_prq;
    __float128 abseps_f2 = abseps / fctr2;
    fctr2 = fctr2 * (8 * constants::pi * state.get_alphaZ());
    Cuba_Int5(data, releps, abseps_f2, res2, err2);

    __float128 res = (res1 * fctr1 + res2 * fctr2) ;
    __float128 err = sqrtq((err1*fctr1)*(err1*fctr1) + (err2*fctr2)*(err2*fctr2)) ;
    print_result(cout, "po1_vertex", res, err, abseps, 6);
}
// void pot1_ver2(quantum_state& state, __float128& abseps, __float128& releps, double& res, double& err)
// {
//     cuba_init5 data5;
//     data5.st = &state;
//     data5.func = &quantum_state::pot1_ver2_I5_prq;
//     __float128 res5, err5;
//     __float128 fctr = 4 * alpha / pow(2 * pi, 7);
//     __float128 fctr_V = 8 * pi * state.get_alphaZ();
//     cuba_init4 data4;
//     data4.st = &state;
//     data4.func = &quantum_state::pot1_ver2_I4_ut;
//     __float128 res4, err4;
//     abseps = abseps / fctr;
//     Cuba_Int5(data5, releps, abseps, res5, err5);
//     Cuba_Int4(data4, releps, abseps, res4, err4);
//     res = (res5 + res4*fctr_V)*fctr *1e6 ;    err = (err5 + err4) * fctr * 1e6 ;
//     // cout << setprecision(16) << (res5 + res4*fctr_V)*fctr *1e6<<
//     //  " ± " << (err5 + err4) * fctr * 1e6  << " x10^6" << "\n" ;
// }

// void pot1(quantum_state& state, double abseps, double releps)
// {
//     __float128res1, err1, res2, err2 ;
//     pot1_ver1(state, abseps, releps, res1, err1);
//     pot1_ver2(state, abseps, releps, res2, err2);
//     cout << "pot1 : " << res1 + res2 <<  " ± " << sqrt(err1*err1 + err2*err2) << " x10^6" << "\n" ;
// }
