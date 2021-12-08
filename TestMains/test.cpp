#include <iostream>
#include "qcd_const.h"
#include "convolution.h"
#include <fstream>
#include <cmath>
using namespace std;


double FF_P(double z, double M){
    double epsilon = std::pow(.3/M,2);
    //double norm = 0.010619;
    return 1./z/std::pow(1-1/z-epsilon/(1-z),2);
}

double FF_NP(double x, double M){
    double norm = 0.0003954134782044029;
    //double norm = 0.010619;
    double a = 25, b = 2.5;
    double a2 = (a+b)*M/4.5-b;
    return std::pow(x,a2-1) * std::pow(1-x,b-1) / norm;
}

double F(double a, double x) {
    return .5* (
        a*a + (-4*a+2*x*x+4*x-6)*std::log(1-x)
      - a*x*(x+2) - x*(x+6) + 4*std::pow(std::log(1-x),2)
    );
}
double FF_NLO(double z, double xmax, double mu2, double M2){
    double r = 1-xmax;
    double D = 1;
    double abar = qcd::CF * qcd::b0/std::log(mu2/qcd::Lambda2);
    double dD = abar*(F(std::log(mu2/M2)-1., xmax)-F(std::log(mu2/M2)-1., 0));
    return (1-dD)/(1-xmax)*(z>xmax) 
          + abar*(1+z*z)/(1.-z)*(std::log(mu2/M2)-2*log(1-z)-1) * (z<xmax);
}

double FF_Peterson(double z, double r){
    return z*pow(1-z,2)/(r*r*z + pow(1-z,2));
}

// https://arxiv.org/pdf/hep-ph/9409316.pdf, HQET Q->PseudoScalar meson
double FF_HF2PS(double z, double r){
    double r2 = r*r, z2 = z*z;
    double r3 = r*r2, z3 = z*z2;
    double r4 = r2*r2, z4 = z2*z2;
    double r5 = r2*r3, z5 = z2*z4;
    double lnr = std::log(r);
    double shape = r*z*pow(1-z,2)/pow(1-(1-r)*z, 6) * (
        6. 
      - 18.*(1.-2*r)*z 
      + (21.-74*r+68*r2)*z2
      - 2*(1.-r)*(6-19*r+18*r2)*z3
      + 3.*pow(1.-r,2) * (1-2*r+2*r2)*z4
    );
    return shape;
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

double FF_Lund(double z, double M){
    //double a = 1.84, b = 0.642;
    double a = 0.89, b = 3.3;
    double Mt2 = M*M+.3+.4;
    return std::pow(1.-z,a)/pow(z,1+b*Mt2) * std::exp(b*Mt2-b*Mt2/z);
}


// https://www.sciencedirect.com/science/article/pii/055032139190597Q



int main() {
    double Q2min = std::pow(1.5+0.33, 2);
    double Q2max = std::pow(10.5/4., 2);
    double tmin = qcd::b0*std::log(std::log(Q2min/qcd::Lambda2));
    double tmax = qcd::b0*std::log(std::log(Q2max/qcd::Lambda2));
    int Nt = 41;
    double dt = (tmax-tmin)/Nt;
    int Nz = 501;
    double zmin = 1e-2;
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
        for (int i=0; i<Nz; i++) dFF[it][i] = 0.;

        if (it.at(0)=='c') {
            double M = qcd::mass_table[it.at(0)];
            for (int i=0; i<Nz; i++) FF[it][i] = 
                z[i]*FF_Lund(z[i], M)*(z[i]>M*M/(M*M+Q2min));
        }
        //else for (int i=0; i<Nz; i++) FF[it][i] = 0.0;
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
        // RK-4
        FFgrids k1, k2, k3, k4, FF1=FF, FF2=FF, FF3=FF;

        Convolve_Valance(t, dlnz, z, zover1mz, FF, dFF, nonsinglets);
        k1 = dFF;        
        for (auto & it : nonsinglets)
            for (int i=0; i<Nz; i++) 
                FF1[it][i] += k1[it][i]*dt/2.;

        Convolve_Valance(t+dt/2., dlnz, z, zover1mz, FF1, dFF, nonsinglets);
        k2 = dFF;
        for (auto & it : nonsinglets)
            for (int i=0; i<Nz; i++) 
                FF2[it][i] += k2[it][i]*dt/2.;

        Convolve_Valance(t+dt/2., dlnz, z, zover1mz, FF2, dFF, nonsinglets);
        k3 = dFF;
        for (auto & it : nonsinglets)
            for (int i=0; i<Nz; i++) 
                FF3[it][i] += k3[it][i]*dt;

        Convolve_Valance(t+dt, dlnz, z, zover1mz, FF3, dFF, nonsinglets);
        k4 = dFF;

        for (auto & it : nonsinglets)
            for (int i=0; i<Nz; i++) 
                FF[it][i] += (k1[it][i] + 2.*k2[it][i] 
                         + 2.*k3[it][i] +    k4[it][i]) * dt/6.;
        
    }

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
                FF[it][i] += (k1[it][i] + 2.*k2[it][i] 
                         + 2.*k3[it][i] +    k4[it][i]) * dt/6.;
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
