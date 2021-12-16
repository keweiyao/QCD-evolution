#ifndef CONVOLUTION_H
#define CONVOLUTION_H

#include <vector>
#include <map>
#include <string>

typedef std::map<std::string, std::vector<double> > FFgrids;

bool ConvolveWithSingular(const double & dlnz, const std::vector<double> & z,
                          const std::vector<double> & zD, const std::vector<double> &R,
                          std::vector<double> & zDxR) ;
bool ConvolveWithSingular2(const double & dlnz, const std::vector<double> & z, const double lnQ2,
                          const std::vector<double> & zD, const std::vector<double> &R,
                          std::vector<double> & zDxR) ;

bool ConvolveWithRegular(const double & dlnz, const std::vector<double> & z,
                          const std::vector<double> & zD, const std::vector<double> &A,
                          std::vector<double> & zDxR);

bool Convolve_Valance(const double & t, const double & dlnz,
                      const std::vector<double> & z, 
                      FFgrids & FF, FFgrids & dFF, std::vector<std::string> flavors, bool med, double jetE);

bool Convolve_Singlets(const double & t, const double & dlnz,
                      const std::vector<double> & z, 
                      FFgrids & FF, FFgrids & dFF, std::vector<std::string> singlets, bool med, double jetE);

#endif
