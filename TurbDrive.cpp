#include <iostream>
#include <cmath>
#include <fftw3.h>
#include <random>
#include <H5Cpp.h>
#include "TurbDrive.h"
#include <iomanip>
#include <sstream>

double compute_k(int i, int j, int l, int n)
{
    int kx = (i <= n/2) ? i : i - n;
    int ky = (j <= n/2) ? j : j - n;
    int kz = (l <= n/2) ? l : l - n;
    return sqrt(kx*kx + ky*ky + kz*kz);
}

std::vector<double> GRF(int n, int k, unsigned int seed)
{
    // Allocate FFTW Complex Array
    fftw_complex* field = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * n * n);
    
    // Initialize random number generator for Gaussian noise
    std::default_random_engine generator(seed);
    std::normal_distribution<double> amplitude_distribution(1.0, 1.0); // Amplitude: Gaussian distribution around unity
    std::uniform_real_distribution<double> phase_distribution(0.0, 2*PI); // Phase: Uniform distribution from 0 to 2pi
    
    // Generate the field in Fourier space
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int l = 0; l < n; ++l) {
                double k_mag = compute_k(i, j, l, n);
                int idx = (i * n + j) * n + l;
                    
                if (k_mag >= k - 1 && k_mag <= k) {
                    double amplitude = std::max(0.0, amplitude_distribution(generator));
                    double phase = phase_distribution(generator);
                    field[idx][0] = amplitude * cos(phase);
                    field[idx][1] = amplitude * sin(phase);
                } else {
                    field[idx][0] = 0.0;
                    field[idx][1] = 0.0;
                }
            }
        }
    }
    
    // Perform inverse FFT to real space
    fftw_complex* real_field = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * n * n);
    fftw_plan plan = fftw_plan_dft_3d(n, n, n, field, real_field, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    
    // Normalize the result
    for (int i = 0; i < n * n * n; ++i) {
        real_field[i][0] /= (n * n * n);
        real_field[i][1] /= (n * n * n);
    }
        
    std::vector<double> Gaussian_Random_Field(n * n * n);
    for (int i = 0; i < n * n * n; ++i) {
        Gaussian_Random_Field[i] = real_field[i][0];  // Store only the real part
    }
    
    
    // Cleanup
    fftw_destroy_plan(plan);
    fftw_free(field);
    fftw_free(real_field);
    
    return Gaussian_Random_Field;
}

void scale_velocity_to_unity(std::vector<double>& pertx, std::vector<double>& perty, std::vector<double>& pertz)
{
    size_t N = pertx.size();
    double current_vrms = 0.0;
    
    // Compute current rms velocity
    for(size_t i = 0; i < N; i++)
    {
        current_vrms += pertx[i] * pertx[i] + perty[i] * perty[i] + pertz[i] * pertz[i];
    }
    current_vrms = std::sqrt(current_vrms / N);
    
    // Compute scale factor
    double scale = 1.0 / current_vrms;
    
    // Multiply scale factor
    for(size_t i = 0; i < N; i++)
    {
        pertx[i] *= scale;
        perty[i] *= scale;
        pertz[i] *= scale;
    }
}

void save_to_hdf5(std::vector<double>& pertx, std::vector<double>& perty, std::vector<double>& pertz, int n, int k, unsigned int seed)
{
    std::stringstream filename;
    filename << "/Users/lvtianning/Documents/SHAO/Turb_Box/zdrv_n" << n << "_k" << k << ".hdf5";
    
    // Create or overwrite 'zdrv.hdf5'
    H5::H5File file(filename.str(), H5F_ACC_TRUNC);
    
    // Define the size of the datasets to (n, n, n)
    hsize_t pert_dims[3] = {static_cast<hsize_t>(n), static_cast<hsize_t>(n), static_cast<hsize_t>(n)};
    H5::DataSpace pert_dataspace(3, pert_dims); // Dataspace for 3D perturbtation array
    
    // Create datasets for pertx, perty, pertz
    H5::DataSet dataset_pertx = file.createDataSet("pertx", H5::PredType::NATIVE_DOUBLE, pert_dataspace);
    H5::DataSet dataset_perty = file.createDataSet("perty", H5::PredType::NATIVE_DOUBLE, pert_dataspace);
    H5::DataSet dataset_pertz = file.createDataSet("pertz", H5::PredType::NATIVE_DOUBLE, pert_dataspace);
    
    // Write the data to the datasets
    dataset_pertx.write(pertx.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_perty.write(perty.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_pertz.write(pertz.data(), H5::PredType::NATIVE_DOUBLE);
    
    // Set attributes
    file.createAttribute("n", H5::PredType::NATIVE_INT, H5::DataSpace(H5S_SCALAR)).write(H5::PredType::NATIVE_INT, &n);
    file.createAttribute("k", H5::PredType::NATIVE_INT, H5::DataSpace(H5S_SCALAR)).write(H5::PredType::NATIVE_INT, &k);
    file.createAttribute("seed", H5::PredType::NATIVE_UINT, H5::DataSpace(H5S_SCALAR)).write(H5::PredType::NATIVE_UINT, &seed);
    
}
