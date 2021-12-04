#include <iostream>
#include "qcd_const.h"
#include "convolution.h"
#include <fstream>
#include <cmath>

int main() {
    double Q2min = std::pow(.4, 2);
    double Q2max = std::pow(100.,2);
    double tmin = qcd::b0*std::log(std::log(Q2min/qcd::Lambda2));
    double tmax = qcd::b0*std::log(std::log(Q2max/qcd::Lambda2));
    int Nt = 201;
    double dt = (tmax-tmin)/Nt;
    int Nz = 201;
    double zmin = 1e-3;
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
            if (it.at(0)=='c' || it.at(0)=='b') FF[it][i] = (1-z[i])*std::pow(z[i], 8)/0.01;
            else if (it.at(0)=='g' || it=="uubar"|| it=="ssbar"|| it=="ddbar") FF[it][i] = 0;
            else FF[it][i] = (1-z[i])*std::pow(z[i], -0.1)/0.6;
            
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
        for (auto & it : nonsinglets){    
            // RK-4
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
