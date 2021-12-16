#include "qcd_const.h"
#include "convolution.h"
#include <cmath>



// convolve zD(z) with P(z) using log-spaced grid
bool ConvolveWithSingular(const double & dlnz, const std::vector<double> & z,
                          const std::vector<double> & zD, const std::vector<double> &R,
                          std::vector<double> & zDxR) {
    zDxR.clear();
    int N = zD.size();
    zDxR.resize(N);
    for (int i=0; i<N; i++) {
        double S1 = 0., S2 = 0.;
        for (int j=i; j<N; j++) S1 += ( zD[N-1-j+i]*R[j] - zD[i]*R[N-1] ) * z[j]/(1.-z[j]);
        for (int j=0; j<i; j++) S2 += z[j]/(1.-z[j]);
        zDxR[i] = (S1 - S2*zD[i]*R[N-1])*dlnz;
    }
    return true;
}

bool ConvolveWithSingular2(const double & dlnz, const std::vector<double> & z, const double lnQ2overLambda2,
                          const std::vector<double> & zD, const std::vector<double> &R,
                          std::vector<double> & zDxR) {
    zDxR.clear();
    int N = zD.size();
    zDxR.resize(N);
    std::vector<double> NL; NL.resize(N);
        for (int i=0; i<N; i++) NL.push_back( 1./(1.+std::max(std::log(1.-z[i]), -.8*lnQ2overLambda2)/lnQ2overLambda2)
                                             );
    for (int i=0; i<N; i++) {
        double S1 = 0., S2 = 0.;
        for (int j=i; j<N; j++) S1 += ( zD[N-1-j+i]*R[j]*NL[j] - zD[i]*R[N-1]*NL[N-1] ) * z[j]/(1.-z[j]);
        for (int j=0; j<i; j++) S2 += z[j]/(1.-z[j])*NL[j];
        zDxR[i] = (S1 - S2*zD[i]*R[N-1])*dlnz;
    }
    return true;
}

bool ConvolveWithRegular(const double & dlnz, const std::vector<double> & z,
                          const std::vector<double> & zD, const std::vector<double> &A,
                          std::vector<double> & zDxR) {
    zDxR.clear();
    int N = zD.size();
    zDxR.resize(N);
    for (int i=0; i<N; i++) {
        double S = 0.;
        for (int j=i; j<N; j++) S += zD[N-1-j+i]*A[j]*z[j];
        zDxR[i] = S*dlnz;
    }
    return true;
}
