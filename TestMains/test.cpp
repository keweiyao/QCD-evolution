#include <iostream>
#include "qcd_const.h"
#include "convolution.h"
#include <fstream>
#include <cmath>
using namespace std;

// https://arxiv.org/pdf/hep-ph/9409316.pdf, HQET Q->PseudoScalar meson
double FF_HF2PS(double z, double r){
    return r*z*pow(1-z,2)/pow(1-(1-r)*z, 6) * (
        6. 
      - 18.*(1.-2*r)*z 
      + (21.-74*r+68*r*r)*pow(z,2)
      - 2*(1.-r)*(6-19*r+18*r*r)*pow(z,3)
      + 3.*pow(1.-r,2) * (1-2*r+2*r*r)*pow(z,4)
    );
}

double FF_HF2VM(double z, double r){
    return 3.*r*z*pow(1-z,2)/pow(1-(1-r)*z, 6) * (
        2. 
      - 2.*(3.-2*r)*z 
      + 3.*(3.-2*r+4*r*r)*pow(z,2)
      - 2*(1.-r)*(4-r+2*r*r)*pow(z,3)
      + pow(1.-r,2) * (3-2*r+2*r*r)*pow(z,4)
    );
}

// https://www.sciencedirect.com/science/article/pii/055032139190597Q
double conter(double x, double a){
    return .5*(a*a+(-4*a+2*x*x+4*x-6)*log(1-x)
          -a*x*(x+2)-x*(x+6) + 4*pow(log(1-x),2)
         );
}
double FF_NLO(double z, double M, double mu2){
    double M2 = M*M;
    double mu = std::sqrt(mu2);
    double R = .3/(M+.3);
    double zmax = 1-M/mu/1.65;//5.2;
    double abar = qcd::CF * qcd::b0/std::log(mu2/qcd::Lambda2);
    double lnQM = std::log(mu2/M2);
    double zmin = M2/(M2+mu2);
    if (z<M2/mu2) return 0;
    return FF_HF2PS(z,R);
    if (z<zmin) return 0.;
    else if (z<zmax) return abar*(1+z*z)/(1.-z)*(lnQM-2*log(1-z)-1);
    else {
        double dx = 1.-zmax;
        double Norm1 = 1.;
        for (int i=0; i<100; i++) {
            double dx = .01*(zmax-zmin);
            double x = zmin+.01*i*(zmax-zmin);
            double aabar = qcd::CF * qcd::b0/std::log(mu2/qcd::Lambda2);
            Norm1 -= dx * aabar*(1+x*x)/(1.-x)*(lnQM-2*log(1-x)-1);
        }
        double Norm2 = 0.;
        for (int i=0; i<20; i++) {
            double x = zmax+.05*dx*i;
            Norm2 += .05*dx * FF_HF2PS(x,R);
        }
        return FF_HF2PS(z,R)*Norm1/Norm2;
    }
}


int main() {
    double Q2min = std::pow(4.9, 2);
    double Q2max = std::pow(92,2);
    double tmin = qcd::b0*std::log(std::log(Q2min/qcd::Lambda2));
    double tmax = qcd::b0*std::log(std::log(Q2max/qcd::Lambda2));
    int Nt = 41;
    double dt = (tmax-tmin)/Nt;
    int Nz = 1001;
    double zmin = .1;
    double zmax = .999;
    double dlnz = std::log(zmax/zmin)/(Nz-1);
    std::vector<double> z, zover1mz; z.clear(), zover1mz.clear();
    for (int i=0; i<Nz; i++) {
        double it = zmin*std::exp(i*dlnz);
        z.push_back(it);
        zover1mz.push_back(it/(1.-it));
    }
    FFgrids FF, dFF;
    std::vector<double> temp; temp.resize(Nz);
    std::vector<std::string> comp{"uv","dv","sv",
                                  "cv","bv",
                                  "g",
                                  "uubar","ddbar","ssbar", 
                                  "ccbar","bbbar"};
    std::vector<std::string> heavy_v{"cv","bv"};
    std::vector<std::string> nonsinglets{"uv","dv","sv",
                                         "cv","bv"};
    std::vector<std::string> singlets{"g",
                                      "uubar","ddbar","ssbar",
                                      "ccbar","bbbar"};
    for (auto & it : comp) {
        FF.insert(std::make_pair(it, temp));
        dFF.insert(std::make_pair(it, temp));
    }
    
    
    for (auto & it : comp) {
        for (int i=0; i<Nz; i++) {
            dFF[it][i] = 0.;
        }
        if (it.at(0)=='b') {
            double m=0.3;
            double M = qcd::mass_table[it.at(0)];
            double r = m/(M+m);
             for (int i=0; i<Nz; i++) FF[it][i] = z[i]*FF_NLO(z[i],M,Q2min);//*( FF_HF2PS(z[i], r, 0) * (z[i]>M*M/Q2min) );
        }
        else {
            for (int i=0; i<Nz; i++) FF[it][i] = 0.0;
        }
    }

    std::ofstream fi("ic.dat");
    for (auto & it : comp) fi << "# " << it << "\t";
    fi << std::endl; 
    for (int i=0; i<Nz; i++) {
        fi << z[i] << "\t";
        for (auto & it : comp) fi << FF[it][i] << "\t";
        fi << std::endl;
    }   

    std::cout << "start decoupled non-singlet evolution" << std::endl;
    for (int i=0; i<Nt; i++) {
        double t = tmin+i*dt;
        double Q2 = qcd::Lambda2*std::exp(std::exp(t/qcd::b0));
        
        for (auto & it : nonsinglets){    
            // RK-4
            double M = qcd::mass_table[it.at(0)];
            FFgrids k1, k2, k3, k4, FF1=FF, FF2=FF, FF3=FF;

            Convolve_Valance(t, dlnz, z, zover1mz, FF, dFF, nonsinglets);
            k1 = dFF;        
            for (int i=0; i<Nz; i++) FF1[it][i] += k1[it][i]*dt/2.;

            Convolve_Valance(t+dt/2., dlnz, z, zover1mz, FF1, dFF, nonsinglets);
            k2 = dFF;
            for (int i=0; i<Nz; i++) FF2[it][i] += k2[it][i]*dt/2.;

            Convolve_Valance(t+dt/2., dlnz, z, zover1mz, FF2, dFF, nonsinglets);
            k3 = dFF;
            for (int i=0; i<Nz; i++) FF3[it][i] += k3[it][i]*dt;

            Convolve_Valance(t+dt, dlnz, z, zover1mz, FF3, dFF, nonsinglets);
            k4 = dFF;

            for (int i=0; i<Nz; i++) FF[it][i] += (k1[it][i] + 2.*k2[it][i] + 2.*k3[it][i] + k4[it][i]) * dt/6.;
            for (int i=0; i<Nz; i++) { if (z[i]<M*M/Q2) FF[it][i]=0.;}
        }
        
    }

    std::ofstream fc("channel.dat");
    std::cout << "start coupled singlets evolution" << std::endl;
    for (int i=0; i<Nt; i++) {
        // RK-4
        double t = tmin+i*dt;
        FFgrids k1, k2, k3, k4, FF1=FF, FF2=FF, FF3=FF;
        Convolve_Singlets(t, dlnz, z, zover1mz, FF, dFF, singlets);
        k1 = dFF;
        for (auto & it : singlets)
            for (int i=0; i<Nz; i++)
                FF1[it][i] += k1[it][i]*dt/2.;

        Convolve_Singlets(t+dt/2., dlnz, z, zover1mz, FF1, dFF, singlets);
        k2 = dFF;
        for (auto & it : singlets)
            for (int i=0; i<Nz; i++)
                FF2[it][i] += k2[it][i]*dt/2.;

        Convolve_Singlets(t+dt/2., dlnz, z, zover1mz, FF2, dFF, singlets);
        k3 = dFF;
        for (auto & it : singlets)
            for (int i=0; i<Nz; i++)
                FF3[it][i] += k3[it][i]*dt;

        Convolve_Singlets(t+dt, dlnz, z, zover1mz, FF3, dFF, singlets);
        k4 = dFF;

        for (auto & it : singlets)
            for (int i=0; i<Nz; i++)
                FF[it][i] += (k1[it][i] + 2.*k2[it][i] + 2.*k3[it][i] + k4[it][i]) * dt/6.;
        fc << qcd::Lambda2*std::exp(std::exp(t/qcd::b0)) << " " << z[30] << " " << FF[singlets[0]][30] << std::endl;
    }


    std::ofstream ff("fs.dat");
    for (auto & it : comp) ff << "# " << it << "\t";
    ff << std::endl; 
    for (int i=0; i<Nz; i++) {
        ff << z[i] << "\t";
        for (auto & it : comp) ff << FF[it][i] << "\t";
        ff << std::endl;
    }   

}
