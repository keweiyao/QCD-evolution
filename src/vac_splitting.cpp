#include <cmath>
#include "qcd_const.h"
#include "vac_splitting.h"
// kt2 * dP/(dkt2 dx)
// Conventions: 
//   P(x) = R(x) / [(1-x)]_+ + A(x) + D(x) delta(1-x)
//   Pij : i->j
double Rqq(double x, double kt2) {
    return qcd::CF * (1.+x*x);
}
double Dqq(double kt2) {
    return qcd::CF_1d5;
}
double Aqg(double x, double kt2) {
    return qcd::CF * (1. + std::pow(1.-x, 2))/x; 
}
double Rgg(double x, double kt2) {
    return qcd::CA_2 * x;
}
double Agg(double x, double kt2) {
    double onemx = 1.-x;
    return qcd::CA_2 * (onemx/x + x*onemx);
}
double Dgg(double kt2, double nf) {
    return qcd::CA_11over6 -  qcd::CA_over9*nf;
}
double Agq(double x, double kt2) {
    return qcd::TR*(pow(x,2) + pow(1.-x,2));
}
