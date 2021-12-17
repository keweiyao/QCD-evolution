#include <fstream>
#include "evolution.h"

double FF_Bowler(double z, double M, double akT){
    //double a = 1.84, b = 0.642;
    //double a = 0.89, b = 3.3;;
    double a = 0.68, b=0.98;
    double Mt2 = M*M+akT*akT;
    //double D = b/z + b*std::log(z) + 1./avgkT2;
    /*double norm = 1;
    for (int i=0; i<100; i++) {
       double dx = (1.-2.*1e-2)/99.;
       double x = 1e-2+dx*i;
       double D = b/x + b*std::log(x) - 1./avgkT2;
       norm += std::pow(1.-x,a)/pow(x,1+b*M2) * std::exp(-b*M2/x) / D * dx;
    }*/
    return std::pow(1.-z,a)/pow(z,1+b*Mt2) * std::exp(b*Mt2*(1.-1./z));
}

double G2Etac(double z){
    double M = qcd::mass_table['c'];
    double Q = 2*M;
    double R0sq = std::pow(0.8,3);
    double a = alphas(Q*Q);
    double norm = R0sq/std::pow(M,3)*a*a/24./M_PI;
    double shape = 3.*z - 2.*z*z + 2.*(1-z)*std::log(1-z);
    return norm * shape;
}

void mdglap::initialize(std::string poi) {
    for (auto & it : comp) for (auto & iit : dZFFdt[it]) iit=0.; // dFF=0
    for (auto & it : comp) for (auto & iit : ZFF[it]) iit=0.; // FF=0
    // for FF
    if (poi.at(0) == 'D')
        for (int iz=0; iz<Nz; iz++)
            ZFF["ccbar"][iz] = z[iz] * FF_Bowler(z[iz], qcd::mass_table['c'], 0.7);
    if (poi.at(0) == 'B')
        for (int iz=0; iz<Nz; iz++)
            ZFF["bbbar"][iz] = z[iz] * FF_Bowler(z[iz], qcd::mass_table['b'], 0.7);
    if (poi.at(0) == 'h') {
        std::vector<std::string> lists({"g","uubar","ddbar","ssbar"});
        for (auto & it : lists){
            std::ifstream fi(std::string("/home/weiyaoke/Documents/New-DGLAP/IC/Qinit/")+it+".dat");
            double zc, a;
            for (int iz=0; iz<Nz; iz++) {
                fi >> zc >> a;
                ZFF[it][iz] = z[iz] *  a;//FF_Bowler(z[iz], 0.0, 0.2);
            }
            fi.close();
        };
    }
    if (poi == std::string("etac") )
        for (int iz=0; iz<Nz; iz++)
            ZFF["g"][iz] = z[iz] * G2Etac(z[iz]);
}

void mdglap::evolve(double tmin, double tmax, double dt, double Nt) {
    /*std::cout << "start decoupled non-singlets evolution ... ";
    for (int i=0; i<Nt; i++) {
	double t = tmin+i*dt;
	// RK-4
	ZFF1=ZFF, ZFF2=ZFF, ZFF3=ZFF;

	Convolve_Valance(t, E, ZFF, dZFFdt);
	k1 = dZFFdt;        
	for (auto & it : nonsinglets)
	    for (int i=0; i<Nz; i++) 
	        ZFF1[it][i] += k1[it][i]*dt/2.;

	Convolve_Valance(t+dt/2., E, ZFF1, dZFFdt);
	k2 = dZFFdt;
	for (auto & it : nonsinglets)
	    for (int i=0; i<Nz; i++) 
	        ZFF2[it][i] += k2[it][i]*dt/2.;

	Convolve_Valance(t+dt/2., E, ZFF2, dZFFdt);
	k3 = dZFFdt;
	for (auto & it : nonsinglets)
	    for (int i=0; i<Nz; i++) 
	        ZFF3[it][i] += k3[it][i]*dt;

	Convolve_Valance(t+dt, E, ZFF3, dZFFdt);
	k4 = dZFFdt;

	for (auto & it : nonsinglets)
	    for (int i=0; i<Nz; i++) 
	        ZFF[it][i] += (k1[it][i] + 2.*k2[it][i] 
	                  + 2.*k3[it][i] +    k4[it][i]) * dt/6.;
    }
    std::cout << "done" << std::endl;*/

    std::cout << "start coupled singlets evolution ... ";
    for (int i=0; i<Nt; i++) {
	// RK-4
	double t = tmin+i*dt;
	ZFF1=ZFF, ZFF2=ZFF, ZFF3=ZFF;

	Convolve_Singlets(t, E, ZFF, dZFFdt);
	k1 = dZFFdt;
	for (auto & it : singlets)
	    for (int i=0; i<Nz; i++)
	        ZFF1[it][i] += k1[it][i]*dt/2.;

	Convolve_Singlets(t+dt/2., E, ZFF1, dZFFdt);
	k2 = dZFFdt;
	for (auto & it : singlets)
	    for (int i=0; i<Nz; i++)
	        ZFF2[it][i] += k2[it][i]*dt/2.;

	Convolve_Singlets(t+dt/2., E, ZFF2, dZFFdt);
	k3 = dZFFdt;
	for (auto & it : singlets)
	    for (int i=0; i<Nz; i++)
	        ZFF3[it][i] += k3[it][i]*dt;

	Convolve_Singlets(t+dt, E, ZFF3, dZFFdt);
	k4 = dZFFdt;

	for (auto & it : singlets)
	    for (int i=0; i<Nz; i++)
	        ZFF[it][i] += (k1[it][i] + 2.*k2[it][i] 
	                  + 2.*k3[it][i] +    k4[it][i]) * dt/6.;
    }
     
    std::cout << "done" << std::endl;
}




void mdglap::Convolve_Valance(double t, double E, FFgrids & FF, FFgrids & dFF) {
    double Q2 = t2Q2(t);
    double asbar = alphas(Q2)/2./M_PI;  // alphas(Q2)/2/pi
    for (auto & it : nonsinglets) {
        // q->q
        double mass = qcd::mass_table[it.at(0)];
        if (mass < 1.0 ) { // treat as light
            std::vector<double> deltaS;
            std::vector<double> Rgrids, Dgrids; Rgrids.clear(); Dgrids.clear();
            for (auto & iz : z) {
                double kt2 = iz*(1-iz)*Q2;
                double MedRqq = mode? MSP.Get("Rqq", {E, kt2, iz}):0.0,
                       MedDqq = mode? MSP.Get("Dqq", {E, Q2}):0.0;
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
                if (iz>xmin) {
                    double MedRQQ = mode? MSP.Get(std::string("R")+(it.at(0))+(it.at(0)), {E, kt2, iz}):0.0;
                    RQQgrids.push_back(RQQ(iz, M2overQ2)+MedRQQ/asbar);
                    AQQgrids.push_back(AQQ(iz, M2overQ2));
                }
                else {
                    RQQgrids.push_back(0.0);
                    AQQgrids.push_back(0.0);
                }
                double MedDQQ = mode? MSP.Get(std::string("D")+(it.at(0))+(it.at(0)), {E, Q2}):0.0;
                DQQgrids.push_back(DQQ(M2overQ2)+MedDQQ/asbar);
            }
            std::vector<double> deltaS, deltaR;
            ConvolveWithSingular(dlnz, z, FF[it], RQQgrids, deltaS);
            ConvolveWithRegular(dlnz, z, FF[it], AQQgrids, deltaR);
            for (int i=0; i<z.size(); i++) dFF[it][i] = deltaS[i] + deltaR[i] + DQQgrids[i]*FF[it][i];
        }
    }
}


void mdglap::Convolve_Singlets(double t, double E, FFgrids & FF, FFgrids & dFF) {
    double Q2 = t2Q2(t);
    double asbar = alphas(Q2)/2./M_PI;  // alphas(Q2)/2/pi
    std::vector<double> RGGgrids, RQQgrids, AGGgrids, AQGgrids, 
                        AGQgrids, DGGgrids, DQQgrids, DGQgrids; 
    RGGgrids.clear(); RQQgrids.clear(); AGGgrids.clear(); AQGgrids.clear(); 
    AGQgrids.clear(); DGGgrids.clear(); DQQgrids.clear(); DGQgrids.clear(); 
    for (auto & iz : z) {
        
        double kt2 = iz*(1-iz)*Q2;
        double MedRgg = mode? MSP.Get("Rgg", {E, kt2, iz}) : 0.0;                 
        RGGgrids.push_back(Rgg(iz, kt2) + MedRgg/asbar);
        double MedAgg = mode? MSP.Get("Agg", {E, kt2, iz}) : 0.0;     
        AGGgrids.push_back(Agg(iz, kt2) + MedAgg/asbar);
        double MedAgq = mode? MSP.Get("Agq", {E, kt2, iz}) : 0.0;  
        AGQgrids.push_back(Agq(iz, kt2) + MedAgq/asbar);
        double MedDgg = mode? MSP.Get("Dgg", {E, Q2}) : 0.0;  
        DGGgrids.push_back(Dgg(kt2) + MedDgg/asbar);
        double MedDgq = mode? MSP.Get("Dgq", {E, Q2}) : 0.0;  
        DGQgrids.push_back(Dgq(kt2) + MedDgq/asbar);
        double MedRqq = mode? MSP.Get("Rqq", {E, kt2, iz}) : 0.0;  
        RQQgrids.push_back(Rqq(iz, kt2) + MedRqq/asbar);
        double MedDqq = mode? MSP.Get("Dqq", {E, Q2}) : 0.0;  
        DQQgrids.push_back(Dqq(kt2) + MedDqq/asbar);
        double MedAqg = mode? MSP.Get("Aqg", {E, kt2, iz}) : 0.0;  
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
                double kt2 = iz*(1-iz)*Q2 - 4.*M2;
                if (kt2>0.) {
                    double MedAGHF = mode? MSP.Get(std::string("Ag")+(iq.at(0)), {E, kt2, iz}) : 0.0;  
                    AGHFgrids.push_back(AgQ(iz, M2overQ2) + MedAGHF/asbar);
                }
                else AGHFgrids.push_back(0.0);
                double MedDGHF = mode? MSP.Get(std::string("Dg")+(iq.at(0)), {E, Q2}) : 0.0;  
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
                double kt2 = iz*(1-iz)*Q2 - std::pow(1-iz,2)*M2;
                if (kt2>0.) {
                    double MedRHH = mode? MSP.Get(std::string("R")+(iq.at(0))+(iq.at(0)), 
                                                 {E, kt2, iz}) : 0.0;
                    RHHgrids.push_back(RQQ(iz, M2overQ2)+MedRHH/asbar);
                    AHHgrids.push_back(AQQ(iz, M2overQ2));
                }
                else {
                    RHHgrids.push_back(0.0);
                    AHHgrids.push_back(0.0);
                }
                double MedDHH = mode? MSP.Get(std::string("D")+(iq.at(0))+(iq.at(0)), 
                                                 {E, Q2}) : 0.0;
                DHHgrids.push_back(DQQ(M2overQ2)+MedDHH/asbar);
            }
            std::vector<double> deltaS, deltaR;
            ConvolveWithSingular(dlnz, z, FF[iq], RHHgrids, deltaS);
            ConvolveWithRegular(dlnz, z, FF[iq], AHHgrids, deltaR);
            for (int i=0; i<z.size(); i++) dFF[iq][i] = deltaS[i] + deltaR[i] + DHHgrids[i]*FF[iq][i];


            // Q->g
            double xmax = Q2/(Q2+M2);
            std::vector<double> AHFGgrids;
            AHFGgrids.clear();
            for (auto & iz : z) {
                double kt2 = iz*(1-iz)*Q2 - std::pow(iz,2)*M2;
                if (kt2>0.) {
                    double MedAHFG = mode? MSP.Get(std::string("A")+iq.at(0)+"g", 
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
}

