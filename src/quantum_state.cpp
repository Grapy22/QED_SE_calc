#include "quantum_state.h"



void quantum_state::ka2jl(const int k , int& j2 , int& l2) 
{   
    if (k){
        j2 = 2 * abs(k) - 1 ;
        l2 = j2 + abs(k) / k ;
    }
    else{
        std::cerr << "Error: k_a != 0 !" << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

quantum_state::quantum_state(int k , int nr , int Z) : k(k) , nr(nr) , Z(Z)
{
    int j2 , l2 ;
    ka2jl(k, j2 , l2) ;
    l = l2 / 2 ;
    j = __float128(j2) / 2 ;
    ka2jl(-k, j2 , l2) ;
    l_ = l2 / 2 ;
    alphaZ = alpha * Z ;
    gamma = sqrtq( k*k - alphaZ * alphaZ ) ; 
    a = alphaZ / sqrtq( powq(nr + gamma ,2) + alphaZ * alphaZ) ;
    en = sqrtq( 1.0Q - a * a  ) ;
    N = 2 * powq(a , 5) / (alphaZ * (alphaZ - k * a) ) ;
    N *= tgammaq(2 * gamma + 1 + nr) / (tgammaq(1 + nr) * powq(tgammaq(2 * gamma + 1), 2)) ;
    C1C2 = C1C2_calc() ;
}
pair <vector<__float128 >, vector<__float128 >> quantum_state::C1C2_calc() const
{
    vector<__float128 > C1(nr+1) , C2(nr+1);
	C1[0] = (-nr) ;
	C2[0] = (alphaZ / a - k);
    __float128  fctr;
	for (int m = 1 ; m < nr +1 ; ++m ){
		fctr = m * (m + 2 * gamma) ;
		C1[m] =  C1[m-1] / fctr * (m - nr);
        C2[m] =  C2[m-1] / fctr * (m - 1 - nr);
	}
	return make_pair(C1 , C2) ;
}

vector<__float128 > quantum_state::Dm_calc (const __float128&  p, const int lmax, const __float128& b) const
{
    const int mmax = nr + lmax;
    const int mmin = -lmax;
    const int m_total = mmax - mmin + 1;
    const int m_start = -mmin;  // Индекс m=0
    const int result_size = nr + 1;

        // Используем только две строки вместо всей матрицы
    vector<__float128 > current_row(m_total, __float128(0));
    vector<__float128 > next_row(m_total, __float128 (0));

        // Предварительные вычисления
    const __float128  a2 = a * a;
    const __float128  a2p2 = a2 + p * p;
    const __float128  fctr = 2 * constants::pi / (a * p);
    const __float128  arctg = atanq(p / a);
    const __float128  inv_p = 1.0Q / p;

    const __float128  ratio_base = (2 * a) / sqrtq(a2p2);

    // Вычисляем начальное значение тяжелой части для самого первого шага (m = mmin)
    const __float128  pwr1_start = mmin + gamma;
    const __float128  pwr2_start = pwr1_start + 1 + b;   
    __float128  running_factor = powq(2 * a, pwr1_start) * tgammaq(pwr2_start) / powq(a2p2, pwr2_start / 2);

        // Заполнение первой строки
    for (int m = mmin; m <= mmax; m++) {
        const int m_idx = m - mmin;
        const __float128  pwr1 = m + gamma;
        const __float128  pwr2 = pwr1 + 1 + b;

        current_row[m_idx] = running_factor * sinq(pwr2 * arctg) * fctr;

        running_factor *= ratio_base * pwr2;
    }

    for (int l__ = 0; l__ < lmax; ++l__) {
    const __float128  l_gamma_b = l__ + 1 + gamma + b;
    const __float128  two_a = 2 * a;

    fill(next_row.begin(), next_row.end(), __float128 (0));

    for (int m = -lmax + 1 + l__; m <= mmax; ++m) {
        const int m_idx = m - mmin;
        if (m_idx >= 1 && m_idx < m_total) {
            next_row[m_idx] = ((l_gamma_b + m) * two_a * current_row[m_idx - 1] - a * current_row[m_idx]) * inv_p;
        }
    }
    // Меняем строки местами
    current_row.swap(next_row);
    }

        // Возвращаем срез
    vector<__float128 > result(result_size);
    for (int i = 0; i <= nr; i++) {
        result[i] = current_row[m_start + i];
    }

    return result;
}

void quantum_state::dirac_pnt_p_dual (const __float128& p, __float128 & g, __float128 & f,
    const __float128&  b , __float128 * dg, __float128 * df ) const
{
    const vector<__float128 >& C1 = C1C2.first;  // Доступ к C1
	const vector<__float128 >& C2 = C1C2.second; // Доступ к C2

    vector<__float128 > D_lk = Dm_calc(p, l, b);
    vector<__float128 > D_l_minus_k = Dm_calc(p, l_, b);

    __float128  sum_g(0) , sum_f(0) ;
	for (int j = 0 ; j<=nr ; ++j){
		sum_g += (C1[j] + C2[j]) * D_lk[j] ;
		sum_f += (C1[j] - C2[j]) * D_l_minus_k[j] ;
	}
	int sign_kappa = (k > 0) - (k < 0); // определение знака
    __float128  factor1 = sqrtq(N * (1.0Q + en)) * p;
    __float128  factor2 = sqrtq(N * (1.0Q - en)) * (-p) * sign_kappa ;
	g = sum_g * factor1 ;
	f = sum_f * factor2 ;

    if (dg or df){
        __float128  dg1 = g * l / p ;
	    __float128  df1 = f * l_ / p ;
        vector<__float128 > D_lk_diff = Dm_calc(p, l + 1, b + 1);
        vector<__float128 > D_l_minus_k_diff = Dm_calc(p, l_ + 1, b + 1);
    	__float128  sum_dg(0) , sum_df(0) ;
		for (int j = 0 ; j<=nr ; j++){
			sum_dg += (C1[j] + C2[j]) * D_lk_diff[j] ;
			sum_df += (C1[j] - C2[j]) * D_l_minus_k_diff[j] ;
		}    
        __float128  dg2 = sum_dg * factor1;
        
	    __float128  df2 = sum_df * factor2;
        *dg = dg1 - dg2 , *df = df1 - df2  ;
    }
}
__float128 quantum_state::Func_prq(pFunc f, const __float128& x, const __float128& y, const __float128& p, const __float128& r,
    const __float128& q) const
{
    __float128 p1 = (p + r) / M_SQRT2q;
    __float128 p2 = (p - r) / M_SQRT2q;
    __float128 ksi = (p1*p1 + p2*p2 - q*q) / (2 * p1 * p2);
    return f(x, y, p1, p2, ksi) * q ;
}