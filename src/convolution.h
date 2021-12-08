#ifndef CONVOLUTION_H
#define CONVOLUTION_H

#include <vector>
#include <map>
#include <string>

typedef std::map<std::string, std::vector<double> > FFgrids;

bool ConvolveWithSingular(const double & dlnz, const double & zover1mz,
                          const std::vector<double> & zD, const std::vector<double> &R,
                          std::vector<double> & zDxR);

bool ConvolveWithRegular(const double & dlnz, const double & z,
                          const std::vector<double> & zD, const std::vector<double> &A,
                          std::vector<double> & zDxR);

bool ConvolveWithDelta(const std::vector<double> & zD, const double D, std::vector<double> & zDxR);

bool Convolve_Valance(const double & t, const double & dlnz,
                      const std::vector<double> & z, const std::vector<double> & zover1mz, 
                      FFgrids & FF, FFgrids & dFF, std::vector<std::string> flavors);
bool Convolve_Singlets(const double & t, const double & dlnz,
                      const std::vector<double> & z, const std::vector<double> & zover1mz, 
                      FFgrids & FF, FFgrids & dFF, std::vector<std::string> singlets);

#endif
