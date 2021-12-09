#include <iostream>
#include "qcd_const.h"
#include "convolution.h"
#include <fstream>
#include <cmath>
using namespace std;


double FF_Lund(double z, double M){
    //double a = 1.84, b = 0.642;
    double a = 0.89, b = 3.3;
    double Mt2 = M*M+.3+.4;
    return std::pow(1.-z,a)/pow(z,1+b*Mt2) * std::exp(b*Mt2-b*Mt2/z);
}


// https://www.sciencedirect.com/science/article/pii/055032139190597Q



int main() {
    double Q2min = std::pow(.4, 2);
    double Q2max = std::pow(6, 2);
    double tmin = qcd::b0*std::log(std::log(Q2min/qcd::Lambda2));
    double tmax = qcd::b0*std::log(std::log(Q2max/qcd::Lambda2));
    int Nt = 51;
    double dt = (tmax-tmin)/Nt;
    int Nz = 4001;
    double zmin = 1e-2;
    double zmax = .9999;
    double dlnz = std::log(zmax/zmin)/(Nz-1);
    std::vector<double> z, zover1mz; z.clear(), zover1mz.clear();
    for (int i=0; i<Nz; i++) {
        double it = zmin*std::exp(i*dlnz);
        z.push_back(it);
        zover1mz.push_back(it/(1.-it));
    }
    
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

    
    for (int i=0; i<2; i++){

    FFgrids FF, dFF;
    for (auto & it : comp) {
        FF.insert(std::make_pair(it, temp));
        dFF.insert(std::make_pair(it, temp));
    }
    bool med = (i==0)? false : true;
    for (auto & it : comp) {
        for (int i=0; i<Nz; i++) dFF[it][i] = 0.;

        if (it == "uv") {
            //std::ifstream ff("/home/weiyaoke/Documents/DGLAP/FF/ic.dat");
            for (int i=0; i<Nz; i++) FF[it][i] = 1.0;
        }
        else for (int i=0; i<Nz; i++) FF[it][i] = 0.0;
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
        std::cout << "Step " << i << ", Q = " << std::sqrt(std::exp(std::exp(t/qcd::b0))*qcd::Lambda2) << std::endl;
        // RK-4
        FFgrids k1, k2, k3, k4, FF1=FF, FF2=FF, FF3=FF;

        Convolve_Valance(t, dlnz, z, zover1mz, FF, dFF, nonsinglets, med);
        k1 = dFF;        
        for (auto & it : nonsinglets)
            for (int i=0; i<Nz; i++) 
                FF1[it][i] += k1[it][i]*dt/2.;

        Convolve_Valance(t+dt/2., dlnz, z, zover1mz, FF1, dFF, nonsinglets, med);
        k2 = dFF;
        for (auto & it : nonsinglets)
            for (int i=0; i<Nz; i++) 
                FF2[it][i] += k2[it][i]*dt/2.;

        Convolve_Valance(t+dt/2., dlnz, z, zover1mz, FF2, dFF, nonsinglets, med);
        k3 = dFF;
        for (auto & it : nonsinglets)
            for (int i=0; i<Nz; i++) 
                FF3[it][i] += k3[it][i]*dt;

        Convolve_Valance(t+dt, dlnz, z, zover1mz, FF3, dFF, nonsinglets, med);
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

        Convolve_Singlets(t, dlnz, z, zover1mz, FF, dFF, singlets, med);
        k1 = dFF;
        for (auto & it : singlets)
            for (int i=0; i<Nz; i++)
                FF1[it][i] += k1[it][i]*dt/2.;

        Convolve_Singlets(t+dt/2., dlnz, z, zover1mz, FF1, dFF, singlets, med);
        k2 = dFF;
        for (auto & it : singlets)
            for (int i=0; i<Nz; i++)
                FF2[it][i] += k2[it][i]*dt/2.;

        Convolve_Singlets(t+dt/2., dlnz, z, zover1mz, FF2, dFF, singlets, med);
        k3 = dFF;
        for (auto & it : singlets)
            for (int i=0; i<Nz; i++)
                FF3[it][i] += k3[it][i]*dt;

        Convolve_Singlets(t+dt, dlnz, z, zover1mz, FF3, dFF, singlets, med);
        k4 = dFF;

        for (auto & it : singlets)
            for (int i=0; i<Nz; i++)
                FF[it][i] += (k1[it][i] + 2.*k2[it][i] 
                         + 2.*k3[it][i] +    k4[it][i]) * dt/6.;
    }


    std::ofstream ff(std::string("./fs-")+std::to_string(i)+std::string(".dat"));
    for (auto & it : comp) ff << "# " << it << "\t";
    ff << std::endl; 
    for (int i=0; i<Nz; i++) {
        ff << z[i] << "\t";
        for (auto & it : comp) ff << FF[it][i] << "\t";
        ff << std::endl;
    }   
    }

}
