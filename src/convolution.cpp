#include "qcd_const.h"
#include "convolution.h"
#include "vac_splitting.h"
#include <cmath>
#include <iostream>

// convolve zD(z) with P(z) using log-spaced grid
bool ConvolveWithSingular(const double & dlnz, const std::vector<double> & zover1mz,
                          const std::vector<double> & zD, const std::vector<double> &R,
                          std::vector<double> & zDxR) {
    zDxR.clear();
    int N = zD.size();
    zDxR.resize(N);
    for (int i=0; i<N; i++) {
        double S1 = 0., S2 = 0.;
        for (int j=i; j<N; j++) S1 += ( zD[N-1-j+i]*R[j] - zD[i]*R[N-1] ) * zover1mz[j];
        for (int j=0; j<i; j++) S2 += zover1mz[j];
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


bool Convolve_Valance(const double & t, const double & dlnz,
                      const std::vector<double> & z, const std::vector<double> & zover1mz, 
                      FFgrids & FF, FFgrids & dFF, std::vector<std::string> flavors) {
    
    double Q2 = qcd::Lambda2*std::exp(std::exp(t/qcd::b0));
    std::vector<double> Rgrids, Dgrids; Rgrids.clear(); Dgrids.clear();
    for (auto & iz : z) {
        double kt2 = iz*(1-iz)*Q2;
        Rgrids.push_back(Rqq(iz, kt2));
        Dgrids.push_back(Dqq(kt2));
    }
    for (auto & it : flavors) {
        std::vector<double> deltaDFF;
        // q->q
        bool status = ConvolveWithSingular(dlnz, zover1mz, FF[it], Rgrids, deltaDFF);
        for (int i=0; i<z.size(); i++) {
            dFF[it][i] = deltaDFF[i] + Dgrids[i]*FF[it][i];
        }
    }
    return true;
}

bool Convolve_Singlets(const double & t, const double & dlnz,
                      const std::vector<double> & z, const std::vector<double> & zover1mz, 
                      FFgrids & FF, FFgrids & dFF, std::vector<std::string> singlets) {
    
    double Q2 = qcd::Lambda2*std::exp(std::exp(t/qcd::b0));
    int Nf = 3;
    std::vector<double> RGGgrids, RQQgrids, AGGgrids, AQGgrids, 
                        AGQgrids, DGGgrids, DQQgrids; 
    RGGgrids.clear(); RQQgrids.clear(); AGGgrids.clear(); AQGgrids.clear(); 
    AGQgrids.clear(); DGGgrids.clear(); DQQgrids.clear();
    for (auto & iz : z) {
        double kt2 = iz*(1-iz)*Q2;
        RGGgrids.push_back(Rgg(iz, kt2));
        AGGgrids.push_back(Agg(iz, kt2));
        AGQgrids.push_back(Agq(iz, kt2));
        AQGgrids.push_back(Aqg(iz, kt2));
        DGGgrids.push_back(Dgg(kt2, Nf));
        RQQgrids.push_back(Rqq(iz, kt2));
        DQQgrids.push_back(Dqq(kt2));
    }
    auto ig = singlets[0]; // gluon
    auto qqbars = singlets;
    qqbars.erase(qqbars.begin()); // qqbar flavors
    // gluon frag
    // g->g
    std::vector<double> deltaS;
    ConvolveWithSingular(dlnz, zover1mz, FF[ig], RGGgrids, deltaS);
    std::vector<double> deltaR;
    ConvolveWithRegular(dlnz, z, FF[ig], AGGgrids, deltaR);
    for (int i=0; i<z.size(); i++)  dFF[ig][i] = deltaS[i] + deltaR[i] + DGGgrids[i]*FF[ig][i];
    // g->q
    for (auto & iq : qqbars) {
        std::vector<double> deltaR;
        ConvolveWithRegular(dlnz, z, FF[iq], AGQgrids, deltaR);
        for (int i=0; i<z.size(); i++)  dFF[ig][i] += deltaR[i];
    }
    // quark frag
    // q->q
    for (auto & iq : qqbars) {
        std::vector<double> deltaS;
        ConvolveWithSingular(dlnz, zover1mz, FF[iq], RQQgrids, deltaS);
        for (int i=0; i<z.size(); i++)  dFF[iq][i] = deltaS[i] + DQQgrids[i]*FF[iq][i];
        // q->g
        std::vector<double> deltaR;
        ConvolveWithRegular(dlnz, z, FF[ig], AQGgrids, deltaR);
        for (int i=0; i<z.size(); i++)  dFF[iq][i] += deltaR[i];
    }

    return true;
}
