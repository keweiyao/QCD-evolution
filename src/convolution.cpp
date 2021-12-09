#include "qcd_const.h"
#include "convolution.h"
#include "vac_splitting.h"
#include "medium_correction.h"
#include "interp_nd.h"
#include <cmath>
#include <iostream>

double Emin = 5., Emax = 1000.;
double kmin = .2, kmax = Emax;
double Zmin = .2/Emax; 
double Zmax = 1.-Zmin;
double E = 20;
MediumCorrections MSP(0,0,0,//10,101,101,
                  std::log(Emin), std::log(Emax), //lnE
                  std::log(kmin), std::log(kmax), //lnkT or lnQ
                  std::log(Zmin/(1.-Zmin)), std::log(Zmax/(1.-Zmax)), //ln[z/(1-z)]
                  "/home/weiyaoke/Documents/SplittingDecomposedGirdCXX/PbPb5020/0-5/g-1.8/"
                 );


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

// convolve zD(z) with P(z) using log-spaced grid
bool ConvolveWithSingular_endpoint(const double & dlnz, const std::vector<double> & z, const double & lnQ2,
                          const std::vector<double> & zD, const std::vector<double> &R,
                          std::vector<double> & zDxR) {
    zDxR.clear();
    int N = zD.size();
    zDxR.resize(N);  
    double Q2 = std::exp(lnQ2);
    std::vector<double> NL; NL.clear();
    for (int j=0; j<N; j++) {
        double kt2 = z[j]*(1.-z[j])*Q2*qcd::Lambda2 - std::pow(1.-z[j], 2)*4.75*4.75;
        if (kt2>2.71828*qcd::Lambda2 && 1.-z[j] > 2.71828/Q2) NL.push_back(1./( 1.+std::log(1-z[j])/lnQ2));
        else if (kt2>2.71828*qcd::Lambda2 && 1.-z[j] < 2.71828/Q2) NL.push_back(1.);
        else NL.push_back(0.);
    }
    for (int i=0; i<N; i++) {
        double S1 = 0., S2 = 0.;
        for (int j=i; j<N; j++) S1 += ( zD[N-1-j+i]*R[j] - zD[i]*R[N-1] ) * z[j]/(1-z[j]) * NL[j];
        for (int j=0; j<i; j++) S2 += z[j]/(1-z[j]) * NL[j];
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
                      const std::vector<double> & z, 
                      FFgrids & FF, FFgrids & dFF, std::vector<std::string> flavors, bool med) {
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
                double MedRqq = med? MSP.Get("Rqq", {E, kt2, iz}):0.0,
                       MedDqq = med? MSP.Get("Dqq", {E, Q2}):0.0;
                Rgrids.push_back((Rqq(iz, kt2)+MedRqq/asbar) );
                Dgrids.push_back(Dqq(kt2) + MedDqq/asbar);
            }
            ConvolveWithSingular(dlnz, z, FF[it], Rgrids, deltaS);
            for (int i=0; i<z.size(); i++) dFF[it][i] = deltaS[i] + Dgrids[i]*FF[it][i];
        } else { // heavy quark
            double M2 = std::pow(mass, 2);
            double M2overQ2 = M2/Q2;
            double xmin = M2/(M2+Q2);
            std::vector<double> RQQgrids, AQQgrids, DQQgrids; 
            RQQgrids.clear(); AQQgrids.clear(); DQQgrids.clear();
            for (auto & iz : z) {
                double kt2 = iz*(1.-iz)*Q2 - std::pow(1.-iz,2)*M2;
                double kt2_over_M2 = kt2/M2;
                if (iz>xmin) {
                    double MedRQQ = med? MSP.Get("Rbb", {E, kt2, iz}):0.0;
                    RQQgrids.push_back(RQQ(iz, M2overQ2)+MedRQQ/asbar);
                    AQQgrids.push_back(AQQ(iz, M2overQ2));
                }
                else {
                    RQQgrids.push_back(0.0);
                    AQQgrids.push_back(0.0);
                }
                double MedDQQ = med? MSP.Get("Dbb", {E, Q2}):0.0;
                DQQgrids.push_back(DQQ(M2overQ2)+MedDQQ/asbar);
            }
            std::vector<double> deltaS, deltaR;
            ConvolveWithSingular_endpoint(dlnz, z, lnQ2, FF[it], RQQgrids, deltaS);
            ConvolveWithRegular(dlnz, z, FF[it], AQQgrids, deltaR);
            for (int i=0; i<z.size(); i++) 
                dFF[it][i] = deltaS[i] + deltaR[i] + DQQgrids[i]*FF[it][i];
        }
    }
    return true;
}

bool Convolve_Singlets(const double & t, const double & dlnz,
                      const std::vector<double> & z, 
                      FFgrids & FF, FFgrids & dFF, std::vector<std::string> singlets, bool med) {
    double lnQ2 = std::exp(t/qcd::b0);
    double Q2 = qcd::Lambda2*std::exp(lnQ2);
    double asbar = qcd::b0/lnQ2;  // alphas(Q2)/2/pi
    std::vector<double> RGGgrids, RQQgrids, AGGgrids, AQGgrids, 
                        AGQgrids, DGGgrids, DQQgrids, DGQgrids; 
    RGGgrids.clear(); RQQgrids.clear(); AGGgrids.clear(); AQGgrids.clear(); 
    AGQgrids.clear(); DGGgrids.clear(); DQQgrids.clear(); DGQgrids.clear(); 
    for (auto & iz : z) {
        double kt2 = iz*(1-iz)*Q2;
        double MedRgg = med? MSP.Get("Rgg", {E, kt2, iz}) : 0.0;                 
        RGGgrids.push_back(Rgg(iz, kt2) + MedRgg/asbar);
        double MedAgg = med? MSP.Get("Agg", {E, kt2, iz}) : 0.0;     
        AGGgrids.push_back(Agg(iz, kt2) + MedAgg/asbar);
        double MedAgq = med? MSP.Get("Agq", {E, kt2, iz}) : 0.0;  
        AGQgrids.push_back(Agq(iz, kt2) + MedAgq/asbar);
        double MedDgg = med? MSP.Get("Dgg", {E, Q2}) : 0.0;  
        DGGgrids.push_back(Dgg(kt2) + MedDgg/asbar);
        double MedDgq = med? MSP.Get("Dgq", {E, Q2}) : 0.0;  
        DGQgrids.push_back(Dgq(kt2) + MedDgq/asbar);
        double MedRqq = med? MSP.Get("Rqq", {E, kt2, iz}) : 0.0;  
        RQQgrids.push_back(Rqq(iz, kt2) + MedRqq/asbar);
        double MedDqq = med? MSP.Get("Dqq", {E, Q2}) : 0.0;  
        DQQgrids.push_back(Dqq(kt2) + MedDqq/asbar);
        double MedAqg = med? MSP.Get("Aqg", {E, kt2, iz}) : 0.0;  
        AQGgrids.push_back(Aqg(iz, kt2) + MedAqg/asbar);
    }
    auto ig = singlets[0]; // gluon
    auto qqbars = singlets;
    qqbars.erase(qqbars.begin()); // qqbar flavors
    // gluon frag
    // g->g
    std::vector<double> deltaS;
    ConvolveWithSingular(dlnz, z, FF[ig], RGGgrids, deltaS);
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
                if (iz<xmin && iz<xmax) {
                    double kt2 = iz*(1-iz)*Q2 - 4*M2;
                    double MedAGHF = med? MSP.Get(std::string("Ag")+std::to_string(iq.at(0)), {E, kt2, iz}) : 0.0;  
                    AGHFgrids.push_back(AgQ(iz, M2overQ2) + MedAGHF/asbar);
                }
                else AGHFgrids.push_back(0.0);
                double MedDGHF = med? MSP.Get(std::string("Dg")+std::to_string(iq.at(0)), {E, Q2}) : 0.0;  
                DGHFgrids.push_back(M2overQ2+ MedDGHF/asbar);
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
            ConvolveWithSingular(dlnz, z, FF[iq], RQQgrids, deltaS);
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
                    double kt2 = iz*(1-iz)*Q2 - std::pow(1-iz,2)*M2;
                    double MedRHH = med? MSP.Get(std::string("R")+std::to_string(iq.at(0))+std::to_string(iq.at(0)), 
                                                 {E, kt2, iz}) : 0.0;
                    RHHgrids.push_back(RQQ(iz, M2overQ2)+MedRHH/asbar);
                    AHHgrids.push_back(AQQ(iz, M2overQ2));
                }
                else {
                    RHHgrids.push_back(0.0);
                    AHHgrids.push_back(0.0);
                }
                double MedDHH = med? MSP.Get(std::string("D")+std::to_string(iq.at(0))+std::to_string(iq.at(0)), 
                                                 {E, 2}) : 0.0;
                DHHgrids.push_back(DQQ(M2overQ2)+MedDHH/asbar);
            }
            std::vector<double> deltaS, deltaR;
            ConvolveWithSingular_endpoint(dlnz, z, lnQ2, FF[iq], RHHgrids, deltaS);
            ConvolveWithRegular(dlnz, z, FF[iq], AHHgrids, deltaR);
            for (int i=0; i<z.size(); i++) dFF[iq][i] = deltaS[i] + deltaR[i] + DHHgrids[i]*FF[iq][i];

            // Q->g
            double xmax = Q2/(Q2+M2);
            std::vector<double> AHFGgrids;
            AHFGgrids.clear();
            for (auto & iz : z) {
                if (iz<xmax) {
                    double kt2 = iz*(1-iz)*Q2 - std::pow(iz,2)*M2;
                    double MedAHFG = med? MSP.Get(std::string("A")+std::to_string(iq.at(0))+"g", 
                                                 {E, kt2, iz}) : 0.0;
                    AHFGgrids.push_back(AQg(iz, M2overQ2)+MedAHFG/asbar);
                }
                else AHFGgrids.push_back(0.0);
            }
            std::vector<double> deltaR2;
            ConvolveWithRegular(dlnz, z, FF[ig], AHFGgrids, deltaR2);
            for (int i=0; i<z.size(); i++)  dFF[iq][i] += deltaR2[i];
        }
    }

    return true;
}
