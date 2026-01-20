// A C++ Code using Mac Low(1999)'s method for turbulence driving, generating the driving hdf5 file in the format ORION2 needed.

#include <iostream>
#include <fftw3.h>
#include <cmath>
#include <random>
#include <vector>
#include <complex>
#include "hdf5.h"
#include "TurbDrive.h"

int main(int argc, const char * argv[]) {
    int n, k;
    unsigned int seed;
    double f_sol;
    // double vrms;
    // double Ein;
    
    // Set grids per dimension n
    std::cout << "Please set grids per dimension n:";
    std::cin >> n;
    
    // Set dimensionless wave number k
    std::cout << "Please set dimensionless wave number k:";
    std::cin >> k;
    
    // Set random seeds
    std::cout << "Please set random seeds for generating pertx, perty and pertz:";
    std::cin >> seed;

    // Choose driving mode
    int mode;
    std::cout << "Select mode (1=single-field Helmholtz mix, 2=two-field normalized mix):";
    std::cin >> mode;

    // Set solenoidal energy fraction
    std::cout << "Please set solenoidal fraction f_sol (0~1, 1=solenoidal, 0=compressive, read README.md to see the differences of the definitions of f_sol between mode 1 and 2):";
    std::cin >> f_sol;
    
    // Set initial rms velocity vrms
    // std::cout << "Please initial rms velocity vrms:";
    // std::cin >> vrms;
    
    // Set energy input Ein
    // std::cout << "Please set energy input Ein:";
    // std::cin >> Ein;
    
    // Generate all the velocity components
    std::vector<double> pertx;
    std::vector<double> perty;
    std::vector<double> pertz;
    GRF(n, k, seed, f_sol, mode, pertx, perty, pertz);
    
    // Scale the velocity to unity
    scale_velocity_to_unity(pertx, perty, pertz);
    
    // Save to HDF5 file
    save_to_hdf5(pertx, perty, pertz, n, k, seed, f_sol, mode);

    return 0;
}
