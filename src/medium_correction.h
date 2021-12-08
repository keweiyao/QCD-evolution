#ifndef MEDIUM_CORRECTION_H
#define MEDIUM_CORRECTION_H

#include "interp_nd.h"
#include <string>
#include <fstream>
#include <map>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>

template <class T>
bool elem_exist(const std::vector<T> & v, T elem){
    return std::find(v.begin(), v.end(), elem) != v.end();
}

typedef std::map<std::string, table_nd> splitting_table_type;

/* Tables: Rqq, Dqq, Aqg,
           Rcc, Dcc, Acg,
           Rbb, Dbb, Abg,
           Rgg, Dgg, Agg,
                     Agq, Agc, Agb;
// 3D tables: R and A[X=lnE, Y=lnkt2, Z=ln(z/(1-z))], 
//            NX=10, NY=101, NZ=101
// e.g. Rqq(D3, {NX, NY, NZ}, {lnXmin, lnYmin, lnZmin}, {lnXmax, lnYmax, lnZmax}),
// 2D tables: D[X=lnE, Y=lnQ2], 
//            NX=10, NY=101
// e.g. Dqq(D2, {NX, NY}, {lnXmin, lnYmin}, {lnXmax, lnYmax}),
// start main body of initilizer
*/

double approx(std::vector<double> x);

class MediumCorrections {
private:
    std::string path_;
    const int D3_, D2_;
    std::vector<std::string> channels_;
    splitting_table_type Tables;

    void prepare_one_table(table_nd & T, std::string filename){
        std::ifstream f(filename.c_str());
        int count = 0;
        double entry;
        std::string line;
        while(getline(f, line)){
            std::istringstream in(line);
            double _;
            if (T.dim()==D3_) in >> _ >> _ >> _ >> entry;
            else if (T.dim()==D2_) in >> _ >> _  >> entry;
            else {
                std::cerr << "Table dim mismatch" << std::endl;
                exit(-1);
            }
            if (count>=T.size()){
                std::cerr << "Table size too large" << std::endl;
                exit(-1);
            }
            T.set_with_linear_index(count, entry);
            count ++;
        }
        if (count<T.size()){
            std::cerr << "Table size too small" << std::endl;
            exit(-1);
        }
    } 
public:
    MediumCorrections(int NX, int NY, int NZ,
                      double lnXmin, double lnXmax,
                      double lnYmin, double lnYmax,
                      double lnZmin, double lnZmax, 
                      std::string path):
    path_(path),
    D3_(3), D2_(2),
    channels_({  "Rqq","Rcc","Rbb","Rgg",
                 "Agg","Agq","Agc","Agb","Aqg","Aqc","Aqb",
                 "Dqq","Dcc","Dbb","Dgg"}) 
    {
        // prepare tables
        for (auto & it : channels_) {
            auto filename = path+"/"+it+".dat";
            std::cout << "Loading table "<< it 
                      << " from " << filename << std::endl;
            if (it.at(0)!='D') {
                table_nd T(D3_, {NX, NY, NZ}, 
                               {lnXmin, lnYmin, lnZmin}, 
                               {lnXmax, lnYmax, lnZmax}, approx);
                prepare_one_table(T, filename);
                Tables.insert({it, T});
            } else {
                table_nd T(D2_, {NX, NY}, 
                               {lnXmin, lnYmin}, 
                               {lnXmax, lnYmax}, approx);
                prepare_one_table(T, filename);
                Tables.insert({it, T});
            }
        }
        std::cout << "done" << std::endl;
    }
    
    double Get(std::string channel, std::vector<double> input){ // e.g. Type=R, a=q, b=q, input=(E, kT2, z)
        // check if this channel exists:
        if (elem_exist<std::string>(channels_, channel)){
            int TableDim = Tables[channel].dim();
            std::vector<double> transformed_input(TableDim);
            transformed_input[0] = std::log(input[0]);
            transformed_input[1] = std::log(input[1]);
            if (TableDim==D3_) transformed_input[2] = std::log(input[2]/(1.-input[2]));
            return Tables[channel].interpolate(transformed_input);
        } else {
            std::cout << "warning: you are looking for a channel that does not exist in medium correction, return 0.0" << std::endl;
            return 0.;
        }
    }
};

#endif 
