#ifndef Evolution
#define Evolution

#include "qcd_const.h"
#include "convolution.h"
#include <string>
#include <cmath>
#include "medium_correction.h"
#include "vac_splitting.h"


class mdglap{
private:
    const int Nz;
    const double zmin, zmax, dlnz;
    bool mode;
    double E;
    std::vector<double> z;
    std::vector<std::string> comp, nonsinglets, singlets;
    FFgrids dZFFdt, ZFF; // main objects
    FFgrids k1, k2, k3, k4, ZFF1, ZFF2, ZFF3; // holder for R-K-4 algorithm
    MediumCorrections MSP;
    void Convolve_Valance(double t, double E, FFgrids & FF, FFgrids & dFF); 
    void Convolve_Singlets(double t, double E, FFgrids & FF, FFgrids & dFF);
public:
    mdglap(int _Nz, double _zmin, double _zmax, std::string medfile):
      Nz(_Nz), zmin(_zmin), zmax(_zmax), dlnz(std::log(zmax/zmin)/(Nz-1)),
      comp({"uv","dv","sv","cv","bv","g","uubar","ddbar","ssbar","ccbar","bbbar"}),
      nonsinglets({"uv","dv","sv","cv","bv"}),
      singlets({"g", "uubar","ddbar","ssbar","ccbar","bbbar"}),
      MSP(10,101,101,
          std::log(5.), std::log(1000.), //lnE
          std::log(.2), std::log(1000.), //lnkT or lnQ
          -std::log(1000./.2-1.), std::log(1000./.2-1.), //ln[z/(1-z)]
          medfile) 
    {
        // assign 
        z.clear();
        for (int i=0; i<Nz; i++) z.push_back(zmin*std::exp(i*dlnz));
        // place hodler for dFF and FF
        std::vector<double> temp; temp.resize(Nz);
        for (auto & it : comp) {
	    ZFF.insert(std::make_pair(it, temp));
	    dZFFdt.insert(std::make_pair(it, temp));
        }
    }
    void initialize(std::string it);
    void setMode(int _mode, double _E) {mode = (_mode==0)? false:true; E=_E; };
    void evolve(double tmin, double tmax, double dt, double Nt);
    void output(std::ofstream & f, std::string description){
        f << "# " <<  description << std::endl; 
        // First line is z: 
	for (int i=0; i<Nz; i++) f << z[i] << " ";
        f << std::endl; 
        // Then go over each channel
	for (auto & it : comp) {
            for (int i=0; i<Nz; i++) f << ZFF[it][i] << " ";
	    f << std::endl;
	}   
    };
    
};

#endif
