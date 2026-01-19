#pragma once

#include <vector>

const double PI = 3.141592653589793;

double compute_k(int i, int j, int l, int n); // Compute magnitude of k

// Create a Gaussian Random Field with n^3 grids and top-hat power spectrum in (k-1)~k
// f_sol sets the solenoidal energy fraction (0=fully compressive, 1=fully solenoidal)
void GRF(int n, int k, unsigned int seed, double f_sol,
		 std::vector<double>& pertx, std::vector<double>& perty, std::vector<double>& pertz);

void scale_velocity_to_unity(std::vector<double>& pertx, std::vector<double>& perty, std::vector<double>& pertz); // Scale the velocity to unity

void save_to_hdf5(std::vector<double>& pertx, std::vector<double>& perty, std::vector<double>& pertz, int n, int k, unsigned int seed, double f_sol); // Save the perturbtation field to the HDF5 file
