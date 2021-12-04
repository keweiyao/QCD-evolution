#include "qcd_const.h"   

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
extern const double qcd::Lambda2 = 0.04;
extern const double qcd::b0 = 2./9.;
extern std::map<char, double> qcd::mass_table{
{'g',0.000},
{'u',0.005},
{'d',0.010},
{'s',0.100},
{'c',1.3},
{'b',4.2},
};


