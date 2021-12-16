#ifndef QCD_CONST_H
#define QCD_CONST_H

#include <map>
#include <string>
namespace qcd{
   extern const double Nc;
   extern const double CA;
   extern const double CF;
   extern const double TR;
   extern const double dA;
   extern const double dF;
   extern const double CA_2;
   extern const double CF_1d5;
   extern const double CA_11over6;
   extern const double CA_over9;
   extern const double Lambda2_q, Lambda2_c, Lambda2_b, NPscale2, alphas_Mz_bar;
   extern const double b0_q, b0_c, b0_b;
   extern std::map<char, double> mass_table;
   extern std::string medium_file;
}

double alphas(double Q2);
double Q22t(double Q2);
double t2Q2(double t);

#endif
