#include <iostream>
#include "qcd_const.h"
#include "convolution.h"
#include <fstream>
#include <cmath>

int main() {
    double Q2min = std::pow(1.,2);
    double Q2max = std::pow(10.,2);
    double tmin = qcd::b0*std::log(std::log(Q2min/qcd::Lambda2));
    double tmax = qcd::b0*std::log(std::log(Q2max/qcd::Lambda2));
    double t = tmin;
    int Nt = 21;
    double dt = (tmax-tmin)/Nt;
    int Nz = 201;
    double zmin = 1e-3;
    double zmax = 1.-zmin;
    double dlnz = std::log(zmax/zmin)/(Nz-1);
    std::vector<double> z, zover1mz; z.clear(), zover1mz.clear();
    for (int i=0; i<Nz; i++) {
        double it = zmin*std::exp(i*dlnz);
        z.push_back(it);
        zover1mz.push_back(it/(1.-it));
    }
    FFgrids FF, dFF;
    std::vector<double> temp; temp.resize(Nz);
    std::vector<std::string> comp{"qv","g","qqbar"};
    std::vector<std::string> nonsinglets{"qv"};
    std::vector<std::string> singlets{"g","qqbar"};
    for (auto & it : comp) {
        FF.insert(std::make_pair(it, temp));
        dFF.insert(std::make_pair(it, temp));
    }
    
    
    for (auto & it : comp) {
        for (int i=0; i<Nz; i++) {
            FF[it][i] = (1-z[i])*std::pow(z[i], -0.1);
            dFF[it][i] = 0.;
        }
    }
    std::ofstream fi("ic.dat");
    for (int i=0; i<Nz; i++) {
        for (auto & it : comp) fi << FF[it][i] << " ";
        fi << std::endl;
    }   

    std::cout << "start non-singlet evolution" << std::endl;
    for (int i=0; i<Nt; i++) {
        std::cout << i << std::endl;
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
                FF[it][i] += (k1[it][i] + 2.*k2[it][i] + 2.*k3[it][i] + k4[it][i]) * dt/6.;
        t+=dt;
    }


    std::cout << "start coupled singlet evolution" << std::endl;
    for (int i=0; i<Nt; i++) {
        std::cout << i << std::endl;
        // RK-4
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
        t+=dt;
    }


    std::ofstream ff("fs.dat");
    for (int i=0; i<Nz; i++) {
        for (auto & it : comp) ff << FF[it][i] << " ";
        ff << std::endl;
    }   

}
