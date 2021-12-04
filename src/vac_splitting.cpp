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
double Dgg(double kt2) {
    return qcd::CA_11over6;
}
double Agq(double x, double kt2) {
    return qcd::TR*(pow(x,2) + pow(1.-x,2));
}
double Dgq(double kt2) {
    return - 1./3.;
}
// Heavy quark related

// (kt2+1mx2m2) * dP/(d(kt2+1mx2m2) dx)
double RQQ(double x, double M2overQ2) {
    return qcd::CF * (1.+x*x) ;
}

double AQQ(double x, double M2overQ2) {
    double onemx = 1.-x;
    return  - qcd::CF * 2. * M2overQ2 ;
}

double DQQ(double M2overQ2) {
    double xmin = 1./(1.+1./M2overQ2);
    double onemx = 1.-xmin;
    double onemx_sq = std::pow(onemx, 2);
    return qcd::CF * (
     onemx + .5*(1.-xmin*xmin) - 2*std::log(1.-xmin)
     + 2. * M2overQ2 * (1.-xmin)
    ); 
}

// (kt2+m2) * dP/(d(kt2+m2) dx)
double AgQ(double x, double M2overQ2) {
    return qcd::TR*(pow(x,2) + pow(1.-x,2) + 2.*M2overQ2); 
}
double DgQ(double M2overQ2) {
    double xmin = .5 * (1. - std::sqrt(1.-4.*M2overQ2));
    double xmax = 1.-xmin;
    return - std::pow(xmax,4)/2. + pow(xmax,3)*2./3. - std::pow(xmax,2)/2. 
           + std::pow(xmin,4)/2. - pow(xmin,3)*2./3. + std::pow(xmin,2)/2. 
           - M2overQ2*(xmax*xmax-xmin*xmin); 
}


double AQg(double x, double M2overQ2) {
    double onemx = 1. - x;
    return qcd::CF * ( (1.+onemx*onemx)/x - 2.*M2overQ2) ;
}
