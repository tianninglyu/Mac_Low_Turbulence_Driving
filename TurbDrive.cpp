#include <iostream>
#include <cmath>
#include <fftw3.h>
#include <random>
#include <H5Cpp.h>
#include "TurbDrive.h"
#include <iomanip>
#include <sstream>
#include <complex>
#include <algorithm>

double compute_k(int i, int j, int l, int n)
{
    int kx = (i <= n/2) ? i : i - n;
    int ky = (j <= n/2) ? j : j - n;
    int kz = (l <= n/2) ? l : l - n;
    return sqrt(kx*kx + ky*ky + kz*kz);
}

void GRF(int n, int k, unsigned int seed, double f_sol,
         std::vector<double>& pertx, std::vector<double>& perty, std::vector<double>& pertz)
{
    const size_t N = static_cast<size_t>(n) * n * n;

    // Clamp solenoidal fraction into [0,1]
    const double fsol = std::clamp(f_sol, 0.0, 1.0);
    const double w_sol = std::sqrt(fsol);
    const double w_comp = std::sqrt(1.0 - fsol);

    // Allocate FFTW complex arrays for three components
    fftw_complex* fieldx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* fieldy = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* fieldz = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* realx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* realy = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* realz = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    // Zero initialize arrays
    for (size_t idx = 0; idx < N; ++idx) {
        fieldx[idx][0] = fieldx[idx][1] = 0.0;
        fieldy[idx][0] = fieldy[idx][1] = 0.0;
        fieldz[idx][0] = fieldz[idx][1] = 0.0;
    }

    // Random generators
    std::default_random_engine generator(seed);
    std::normal_distribution<double> amplitude_distribution(1.0, 1.0); // Amplitude around unity
    std::uniform_real_distribution<double> phase_distribution(0.0, 2 * PI); // Phase uniform in [0, 2pi]

    // Build Fourier-space vector field with controlled solenoidal/compressive mix
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int l = 0; l < n; ++l) {
                const int idx = (i * n + j) * n + l;

                // Compute k components
                const int kx = (i <= n / 2) ? i : i - n;
                const int ky = (j <= n / 2) ? j : j - n;
                const int kz = (l <= n / 2) ? l : l - n;
                const double k_mag = std::sqrt(static_cast<double>(kx * kx + ky * ky + kz * kz));

                // Skip k=0 and modes outside the top-hat shell
                if (k_mag < 1e-12 || k_mag < k - 1 || k_mag > k) {
                    continue;
                }

                // Random complex vector amplitude
                const double amp_x = std::max(0.0, amplitude_distribution(generator));
                const double amp_y = std::max(0.0, amplitude_distribution(generator));
                const double amp_z = std::max(0.0, amplitude_distribution(generator));
                const double phase_x = phase_distribution(generator);
                const double phase_y = phase_distribution(generator);
                const double phase_z = phase_distribution(generator);

                const std::complex<double> ax = amp_x * std::complex<double>(std::cos(phase_x), std::sin(phase_x));
                const std::complex<double> ay = amp_y * std::complex<double>(std::cos(phase_y), std::sin(phase_y));
                const std::complex<double> az = amp_z * std::complex<double>(std::cos(phase_z), std::sin(phase_z));

                // Projection onto k-hat
                const double hx = kx / k_mag;
                const double hy = ky / k_mag;
                const double hz = kz / k_mag;
                const std::complex<double> dot = ax * hx + ay * hy + az * hz;

                const std::complex<double> a_para_x = dot * hx;
                const std::complex<double> a_para_y = dot * hy;
                const std::complex<double> a_para_z = dot * hz;

                const std::complex<double> a_perp_x = ax - a_para_x;
                const std::complex<double> a_perp_y = ay - a_para_y;
                const std::complex<double> a_perp_z = az - a_para_z;

                // Blend solenoidal and compressive parts
                const std::complex<double> a_mix_x = w_sol * a_perp_x + w_comp * a_para_x;
                const std::complex<double> a_mix_y = w_sol * a_perp_y + w_comp * a_para_y;
                const std::complex<double> a_mix_z = w_sol * a_perp_z + w_comp * a_para_z;

                fieldx[idx][0] = a_mix_x.real();
                fieldx[idx][1] = a_mix_x.imag();
                fieldy[idx][0] = a_mix_y.real();
                fieldy[idx][1] = a_mix_y.imag();
                fieldz[idx][0] = a_mix_z.real();
                fieldz[idx][1] = a_mix_z.imag();
            }
        }
    }

    // Inverse FFT to real space
    fftw_plan plan_x = fftw_plan_dft_3d(n, n, n, fieldx, realx, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan plan_y = fftw_plan_dft_3d(n, n, n, fieldy, realy, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan plan_z = fftw_plan_dft_3d(n, n, n, fieldz, realz, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(plan_x);
    fftw_execute(plan_y);
    fftw_execute(plan_z);

    // Normalize and copy real parts
    pertx.resize(N);
    perty.resize(N);
    pertz.resize(N);
    for (size_t idx = 0; idx < N; ++idx) {
        pertx[idx] = realx[idx][0] / static_cast<double>(N);
        perty[idx] = realy[idx][0] / static_cast<double>(N);
        pertz[idx] = realz[idx][0] / static_cast<double>(N);
    }

    // Cleanup
    fftw_destroy_plan(plan_x);
    fftw_destroy_plan(plan_y);
    fftw_destroy_plan(plan_z);
    fftw_free(fieldx);
    fftw_free(fieldy);
    fftw_free(fieldz);
    fftw_free(realx);
    fftw_free(realy);
    fftw_free(realz);
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

void save_to_hdf5(std::vector<double>& pertx, std::vector<double>& perty, std::vector<double>& pertz, int n, int k, unsigned int seed, double f_sol)
{
    std::stringstream filename;
    filename << "zdrv_n" << n << "_k" << k << ".hdf5";
    
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
    file.createAttribute("f_sol", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(H5S_SCALAR)).write(H5::PredType::NATIVE_DOUBLE, &f_sol);
    
}
