#include "qcd_const.h"
#include "convolution.h"
#include "vac_splitting.h"
#include <cmath>
#include <iostream>

double Emin = 5., Emax = 1000.;
double kmin = .2, kmax = Emaxs;
double Zmin = .2/Emax; double Zmax = 1.-Zmin;

MediumCorrections MSP(10,101,101,
                  std::log(Emin), std::log(Emax), //lnE
                  std::log(kmin*kmin), std::log(kmax*kmax), //lnkT2 or lnQ2
                  std::log(Zmin/(1.-Zmin)), std::log(Zmax/(1.-Zmax)), //ln[z/(1-z)]
                  "./Tables/"
                 );

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

// convolve zD(z) with P(z) using log-spaced grid
bool ConvolveWithSingular_endpoint(const double & dlnz, const std::vector<double> & zover1mz, const double & lnQ2,
                          const std::vector<double> & zD, const std::vector<double> &R,
                          std::vector<double> & zDxR) {
    zDxR.clear();
    int N = zD.size();
    zDxR.resize(N);  
    std::vector<double> NL; NL.clear();
    for (int j=0; j<N; j++) {
        double z = zover1mz[j]/(1.+zover1mz[j]);
        NL.push_back(1./(
             1.+std::log(std::max(1.-z, 2.71828/std::exp(lnQ2))
                        )/lnQ2
                      )
                    );
    }
    for (int i=0; i<N; i++) {
        double S1 = 0., S2 = 0.;
        for (int j=i; j<N; j++) S1 += ( zD[N-1-j+i]*R[j] - zD[i]*R[N-1] ) * zover1mz[j] * NL[j];
        for (int j=0; j<i; j++) S2 += zover1mz[j] * NL[j];
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
    double lnQ2 = std::exp(t/qcd::b0);
    double Q2 = qcd::Lambda2*std::exp(lnQ2);
    double asbar = qcd::b0/lnQ2;  // alphas(Q2)/2/pi
    for (auto & it : flavors) {
        
        // q->q
        double mass = qcd::mass_table[it.at(0)];
        if (mass < 1.0 ) { // treat as light
            std::vector<double> deltaS;
            std::vector<double> Rgrids, Dgrids; Rgrids.clear(); Dgrids.clear();
            for (auto & iz : z) {
                double kt2 = iz*(1-iz)*Q2;
                double MedRqq = MSP.Get("Rqq", {E, kt2, iz}),
                       MedDqq = MSP.Get("Dqq", {E, kt2});
                Rgrids.push_back(Rqq(iz, kt2) + MedRqq/asbar );
                Dgrids.push_back(Dqq(kt2) + MedDqq/asbar);
            }
            ConvolveWithSingular_endpoint(dlnz, zover1mz, lnQ2, FF[it], Rgrids, deltaS);
            for (int i=0; i<z.size(); i++) dFF[it][i] = deltaS[i] + Dgrids[i]*FF[it][i];
        } else { // heavy quark
            double M2 = std::pow(mass, 2);
            double M2overQ2 = M2/Q2;
            double xmin = M2/(M2+Q2);
            std::vector<double> RQQgrids, AQQgrids, DQQgrids; 
            RQQgrids.clear(); AQQgrids.clear(); DQQgrids.clear();
            for (auto & iz : z) {
                double kt2_over_M2 = iz*(1.-iz)*Q2/M2 - std::pow(1.-iz,2);
                if (iz>xmin) {
                    double MedRQQ = MSP.Get("Rqq", {E, kt2, iz}),
                           MedDQQ = MSP.Get("Dqq", {E, kt2});
                    RQQgrids.push_back(RQQ(iz, M2overQ2) + MedRQQ/asbar );
                    AQQgrids.push_back(AQQ(iz, M2overQ2) + MedDQQ/asbar);
                }
                else {
                    RQQgrids.push_back(0.0);
                    AQQgrids.push_back(0.0);
                }
                DQQgrids.push_back(DQQ(M2overQ2));
            }
            std::vector<double> deltaS, deltaR;
            ConvolveWithSingular_endpoint(dlnz, zover1mz, lnQ2, FF[it], RQQgrids, deltaS);
            ConvolveWithRegular(dlnz, z, FF[it], AQQgrids, deltaR);
            for (int i=0; i<z.size(); i++) 
                dFF[it][i] = deltaS[i] + deltaR[i] + DQQgrids[i]*FF[it][i];
        }
    }
    return true;
}

bool Convolve_Singlets(const double & t, const double & dlnz,
                      const std::vector<double> & z, const std::vector<double> & zover1mz, 
                      FFgrids & FF, FFgrids & dFF, std::vector<std::string> singlets) {
    double lnQ2 = std::exp(t/qcd::b0);
    double Q2 = qcd::Lambda2*std::exp(lnQ2);
    std::vector<double> RGGgrids, RQQgrids, AGGgrids, AQGgrids, 
                        AGQgrids, DGGgrids, DQQgrids, DGQgrids; 
    RGGgrids.clear(); RQQgrids.clear(); AGGgrids.clear(); AQGgrids.clear(); 
    AGQgrids.clear(); DGGgrids.clear(); DQQgrids.clear(); DGQgrids.clear(); 
    for (auto & iz : z) {
        double kt2 = iz*(1-iz)*Q2;
        RGGgrids.push_back(Rgg(iz, kt2));
        AGGgrids.push_back(Agg(iz, kt2));
        AGQgrids.push_back(Agq(iz, kt2));
        AQGgrids.push_back(Aqg(iz, kt2));
        DGGgrids.push_back(Dgg(kt2));
        RQQgrids.push_back(Rqq(iz, kt2));
        DQQgrids.push_back(Dqq(kt2));
        DGQgrids.push_back(Dgq(kt2));
    }
    auto ig = singlets[0]; // gluon
    auto qqbars = singlets;
    qqbars.erase(qqbars.begin()); // qqbar flavors
    // gluon frag
    // g->g
    std::vector<double> deltaS;
    ConvolveWithSingular_endpoint(dlnz, zover1mz, lnQ2, FF[ig], RGGgrids, deltaS);
    std::vector<double> deltaR;
    ConvolveWithRegular(dlnz, z, FF[ig], AGGgrids, deltaR);
    for (int i=0; i<z.size(); i++)  dFF[ig][i] = deltaS[i] + deltaR[i] + DGGgrids[i]*FF[ig][i];
    // g->q
    for (auto & iq : qqbars) {
        double mass = qcd::mass_table[iq.at(0)];
        if (mass < 1.0 ) { // treat as light
            std::vector<double> deltaR;
            ConvolveWithRegular(dlnz, z, FF[iq], AGQgrids, deltaR);
            for (int i=0; i<z.size(); i++)  dFF[ig][i] += deltaR[i] + DGQgrids[i]*FF[ig][i];;
        } else { // g to Q + Qbar
            double M2 = mass*mass;
            if (Q2<4.*M2) continue; // below threshold
            double M2overQ2 = M2/Q2;
            double xmin = .5 * (1. - std::sqrt(1.-4.*M2overQ2));
            double xmax = 1.-xmin;
            std::vector<double> AGHFgrids, DGHFgrids; 
            AGHFgrids.clear(); DGHFgrids.clear();
            for (auto & iz : z) {
                if (iz<xmin && iz<xmax) AGHFgrids.push_back(AgQ(iz, M2overQ2));
                else AGHFgrids.push_back(0.0);
                DGHFgrids.push_back(M2overQ2);
            }
            std::vector<double> deltaR;
            ConvolveWithRegular(dlnz, z, FF[iq], AGHFgrids, deltaR);
            for (int i=0; i<z.size(); i++)  dFF[ig][i] += deltaR[i] + DGHFgrids[i]*FF[ig][i];
        }
    }
    // quark frag
    for (auto & iq : qqbars) {
        double mass = qcd::mass_table[iq.at(0)];
        if (mass < 1.0 ) { // treat as light
            // q->q
            std::vector<double> deltaS;
            ConvolveWithSingular_endpoint(dlnz, zover1mz, lnQ2, FF[iq], RQQgrids, deltaS);
            for (int i=0; i<z.size(); i++)  dFF[iq][i] = deltaS[i] + DQQgrids[i]*FF[iq][i];
            // q->g
            std::vector<double> deltaR;
            ConvolveWithRegular(dlnz, z, FF[ig], AQGgrids, deltaR);
            for (int i=0; i<z.size(); i++)  dFF[iq][i] += deltaR[i];
        } else { 
            // Q->Q
            double M2 = std::pow(mass, 2);
            double M2overQ2 = M2/Q2;
            double xmin = M2/(M2+Q2);
            std::vector<double> RHHgrids, AHHgrids, DHHgrids; 
            RHHgrids.clear(); AHHgrids.clear(); DHHgrids.clear();
            for (auto & iz : z) {
                if (iz>xmin) {
                    RHHgrids.push_back(RQQ(iz, M2overQ2));
                    AHHgrids.push_back(AQQ(iz, M2overQ2));
                }
                else {
                    RHHgrids.push_back(0.0);
                    AHHgrids.push_back(0.0);
                }
                DHHgrids.push_back(DQQ(M2overQ2));
            }
            std::vector<double> deltaS, deltaR;
            ConvolveWithSingular_endpoint(dlnz, zover1mz, lnQ2, FF[iq], RHHgrids, deltaS);
            ConvolveWithRegular(dlnz, z, FF[iq], AHHgrids, deltaR);
            for (int i=0; i<z.size(); i++) dFF[iq][i] = deltaS[i] + deltaR[i] + DHHgrids[i]*FF[iq][i];

            // Q->g
            double xmax = Q2/(Q2+M2);
            std::vector<double> AHFGgrids;
            AHFGgrids.clear();
            for (auto & iz : z) {
                if (iz<xmax) AHFGgrids.push_back(AQg(iz, M2overQ2));
                else AHFGgrids.push_back(0.0);
            }
            std::vector<double> deltaR2;
            ConvolveWithRegular(dlnz, z, FF[ig], AHFGgrids, deltaR2);
            for (int i=0; i<z.size(); i++)  dFF[iq][i] += deltaR2[i];
        }
    }

    return true;
}
