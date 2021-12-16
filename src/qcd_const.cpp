#include "qcd_const.h"   
#include <cmath>
extern const double qcd::Nc = 3.;
extern const double qcd::CA = qcd::Nc;
extern const double qcd::CF = (qcd::Nc*qcd::Nc - 1)/(2.*qcd::Nc);
extern const double qcd::TR = .5;
extern const double qcd::dA = qcd::Nc*qcd::Nc - 1.;
extern const double qcd::dF = qcd::Nc;
extern const double qcd::CA_2 = qcd::CA * 2.;
extern const double qcd::CF_1d5 = qcd::CF * 1.5;
extern const double qcd::CA_11over6 = qcd::CA * 11./6.;
extern const double qcd::CA_over9 = qcd::CA / 9.;

extern std::map<char, double> qcd::mass_table{
{'g',0.000},
{'u',0.005},
{'d',0.010},
{'s',0.100},
{'c',1.3},
{'b',4.5},
};
double Mc = qcd::mass_table['c'],
       Mb = qcd::mass_table['b'],
       Mz = 91.7;  

extern const double qcd::b0_q = 2./(11.-2./3.*3);
extern const double qcd::b0_c = 2./(11.-2./3.*4);
extern const double qcd::b0_b = 2./(11.-2./3.*5);
extern const double qcd::alphas_Mz_bar = 0.118 / 2. / M_PI;
extern const double qcd::Lambda2_b = std::pow(Mz,2)/std::exp(qcd::b0_b/qcd::alphas_Mz_bar);
double a5_Mb_bar = qcd::b0_b/std::log(std::pow(2*Mb,2)/qcd::Lambda2_b);
extern const double qcd::Lambda2_c = std::pow(2*Mb,2)/std::exp(qcd::b0_c/a5_Mb_bar);
double a4_Mc_bar = qcd::b0_c/std::log(std::pow(2*Mc,2)/qcd::Lambda2_c);
extern const double qcd::Lambda2_q = std::pow(2*Mc,2)/std::exp(qcd::b0_q/a4_Mc_bar);
extern const double qcd::NPscale2 = qcd::Lambda2_q * 4.0;

const double tc = qcd::b0_q*std::log(std::log(4*Mc*Mc/qcd::Lambda2_q));
const double tb = qcd::b0_c*std::log(std::log(4*Mb*Mb/qcd::Lambda2_c))
                - qcd::b0_c*std::log(std::log(4*Mc*Mc/qcd::Lambda2_c))
                + qcd::b0_q*std::log(std::log(4*Mc*Mc/qcd::Lambda2_q));

double alphas(double Q2) {
    Q2 = std::max(qcd::NPscale2, Q2);
    if (Q2<4*Mc*Mc) return 2.*M_PI*qcd::b0_q/std::log(Q2/qcd::Lambda2_q);
    else if (Q2<4*Mb*Mb) return 2.*M_PI*qcd::b0_c/std::log(Q2/qcd::Lambda2_c);
    else return 2.*M_PI*qcd::b0_b/std::log(Q2/qcd::Lambda2_b);
}

double Q22t(double Q2) {
    Q2 = std::max(qcd::NPscale2, Q2);
    if (Q2<4*Mc*Mc) 
        return qcd::b0_q*std::log(std::log(Q2/qcd::Lambda2_q));
    else if (Q2<4*Mb*Mb) 
        return qcd::b0_c*std::log(std::log(Q2/qcd::Lambda2_c))
             - qcd::b0_c*std::log(std::log(4*Mc*Mc/qcd::Lambda2_c))
             + qcd::b0_q*std::log(std::log(4*Mc*Mc/qcd::Lambda2_q));
    else 
        return qcd::b0_b*std::log(std::log(Q2/qcd::Lambda2_b))
             - qcd::b0_b*std::log(std::log(4*Mb*Mb/qcd::Lambda2_b))
             + qcd::b0_c*std::log(std::log(4*Mb*Mb/qcd::Lambda2_c))
             - qcd::b0_c*std::log(std::log(4*Mc*Mc/qcd::Lambda2_c))
             + qcd::b0_q*std::log(std::log(4*Mc*Mc/qcd::Lambda2_q)); 

}

double t2Q2(double t) {
    if (t<tc) 
        return qcd::Lambda2_q*std::exp(std::exp(t/qcd::b0_q));
    else if (t<tb) {
        t = t-tc + qcd::b0_c*std::log(std::log(4*Mc*Mc/qcd::Lambda2_c));
        return qcd::Lambda2_c*std::exp(std::exp(t/qcd::b0_c));
    }
    else {
        t = t-tb + qcd::b0_b*std::log(std::log(4*Mb*Mb/qcd::Lambda2_b));
        return qcd::Lambda2_b*std::exp(std::exp(t/qcd::b0_b));
    }
}
