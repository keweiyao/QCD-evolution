#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include "qcd_const.h"
#include "evolution.h"
#include <fstream>
#include <cmath>
#include <experimental/filesystem>


std::string D2S(double value) {
    return "s";
    std::stringstream sstr;
    sstr << std::setprecision(0) << value;
    return sstr.str();
}


int main(int argc, char * argv[]) {
    std::vector<std::string> POIs{"B"}; // particle of interests, inclusive hadron, D meson, B meson

    //std::cout << "Alphas check: " << alphas(1.0) << " "  << alphas(10*10) << " " << alphas(91.7*91.7) << std::endl;
    std::string A = argv[1], B = argv[2], Sqrts = argv[3], cen=argv[4], gs=argv[5];
    double sqrts = atof(argv[3]);
    std::string OutputFolder(std::string(argv[6])+"/"+A+B+Sqrts+"/"+cen+"/g-"+gs+"/");
    std::string TableFolder(std::string(argv[7])+"/"+A+B+Sqrts+"/"+cen+"/g-"+gs+"/");
    std::cout << "System: " << A << " + " << B << " at " << sqrts << std::setprecision(1) << "GeV" << std::endl;
    std::cout << "Loading med-splitting table from " << TableFolder << std::endl;
    std::cout << "Output to " << OutputFolder << std::endl;

    // Create output folder
    std::experimental::filesystem::path OutputPath = std::experimental::filesystem::path(OutputFolder.c_str());
    //if (!std::experimental::filesystem::exists(OutputPath)) std::experimental::filesystem::create_directory(OutputPath.c_str());

    int NE = 15;
    double Emin = 5.0, Emax = 500., ktmin = 0.2;
    std::vector<double> Es; Es.clear();
    std::cout << "Modified fragmentation function will be evaluated at: E [GeV] = " << std::endl;
    for (int i=0; i<NE; i++) {
        Es.push_back(Emin * std::exp(i*std::log(Emax/Emin)/(NE-1)));
        std::cout << Es.back() << " ";
    }
    std::cout << std::endl;

    int Nz = 3001;
    double zmin=1e-2, zmax=.9999;
    mdglap solver(Nz, zmin, zmax, TableFolder);

    // Evolve the IC to starting scale
    /*for (auto & poi : POIs) {
        std::stringstream sf1, sf2, sf3;
        sf1  << "./" << poi << "-Q0.dat"; // initial condition
        sf2  << "./" << poi << "-Q-vac.dat"; // vac frag
        std::ofstream f1(sf1.str()), f2(sf2.str());
            double tmin = Q22t(0.16);
            double tmax = Q22t(1.3);
            int Nt = 50;//int(std::abs(tmax-tmin)/a_reasonable_dt);
            double dt = (tmax-tmin)/Nt;
            solver.initialize(poi);
            solver.setMode(0, 0); // 0 vacuum and 1 for vac+med
            solver.evolve(tmin, tmax, dt, Nt);
            solver.output(f2, "IC");
        f1.close(); f2.close();

    }

    


    return 0; */
    // a reasonable choice of evolution 
    // consider the typical medium scale is 2.0 GeV
    // then dt = t[2.0*1.1]-t[2.0]
    double TypicalMedScale = 2.0;
    double a_reasonable_dt = Q22t(std::pow(1.1*TypicalMedScale,2)) - Q22t(std::pow(TypicalMedScale,2));
    for (auto & poi : POIs) {
        std::stringstream sf1, sf2, sf3;
        sf1 << OutputPath.c_str() << "/" << poi << "-Q0.dat"; // initial condition
        sf2 << OutputPath.c_str() << "/" << poi << "-Q-vac.dat"; // vac frag
        sf3 << OutputPath.c_str() << "/" << poi << "-Q-vac+med.dat"; // vac+med frag
        std::ofstream f2(sf2.str()), f3(sf3.str());

        double mass_scale = 0.0;
        if (poi.at(0)=='D') mass_scale = qcd::mass_table['c'];
        if (poi.at(0)=='B') mass_scale = qcd::mass_table['b'];
	if (cen==std::string("0-5") && gs==std::string("1.6")){ 
			std::cout << "Run vac" << std::endl;
           for (int ie=0; ie<Es.size(); ie++) {
            double jetE = Es[ie];
            double Q2min = std::pow(2*ktmin, 2);
            double Q2max = (std::pow(jetE, 2) + std::pow(mass_scale,2));
            double tmin = Q22t(Q2min);
            double tmax = Q22t(Q2max);
            int Nt = int(std::abs(tmax-tmin)/a_reasonable_dt);
            double dt = (tmax-tmin)/Nt;
            std::cout << poi << ", E[GeV] = " << jetE <<", Tsteps = " << Nt << std::endl;   
            // Vac evolution
            solver.initialize(poi);
            solver.setMode(0, jetE); // 0 vacuum and 1 for vac+med
            solver.evolve(tmin, tmax, dt, Nt);
            solver.output(f2, D2S(jetE));
          }
	}
        for (int ie=0; ie<Es.size(); ie++) {
            double jetE = Es[ie];
            double Q2min = std::pow(2*ktmin, 2);
            double Q2max = (std::pow(jetE, 2) + std::pow(mass_scale,2));
            double tmin = Q22t(Q2min);
            double tmax = Q22t(Q2max);
            int Nt = int(std::abs(tmax-tmin)/a_reasonable_dt);
            double dt = (tmax-tmin)/Nt;
            std::cout << poi << ", E[GeV] = " << jetE <<", Tsteps = " << Nt << std::endl;   
            // Vac+med evolution
            solver.initialize(poi);
            solver.setMode(1, jetE); // 0 vacuum and 1 for vac+med
            solver.evolve(tmin, tmax, dt, Nt);
            solver.output(f3, D2S(jetE));	    
        }  
        f2.close(); f3.close();
    }

    return 0;
}
