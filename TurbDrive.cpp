#include <iostream>
#include <cmath>
#include <fftw3.h>
#include <random>
#include "hdf5.h"
#include "TurbDrive.h"
#include <iomanip>
#include <sstream>
#include <complex>
#include <algorithm>

void GRF(int n, int k, unsigned int seed, double f_sol, int mode,
         std::vector<double>& pertx, std::vector<double>& perty, std::vector<double>& pertz)
{
    const size_t N = static_cast<size_t>(n) * n * n;

    const double fsol = std::clamp(f_sol, 0.0, 1.0);

    double w_sol_mix = 0.0;
    double w_comp_mix = 0.0;
    if (mode == 1) {
        w_sol_mix = fsol;
        w_comp_mix = 1.0 - fsol;
    } else {
        w_sol_mix = std::sqrt(fsol);
        w_comp_mix = std::sqrt(1.0 - fsol);
    }
    
    auto build_field = [&](double w_sol, double w_comp, unsigned int field_seed,
                          std::vector<double>& outx, std::vector<double>& outy, std::vector<double>& outz)
    {
        fftw_complex* fieldx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        fftw_complex* fieldy = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        fftw_complex* fieldz = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        fftw_complex* realx  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        fftw_complex* realy  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        fftw_complex* realz  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

        for (size_t idx = 0; idx < N; ++idx) {
            fieldx[idx][0] = fieldx[idx][1] = 0.0;
            fieldy[idx][0] = fieldy[idx][1] = 0.0;
            fieldz[idx][0] = fieldz[idx][1] = 0.0;
        }

        std::mt19937 generator(field_seed);
        std::normal_distribution<double> amplitude_distribution(1.0, 1.0);
        std::uniform_real_distribution<double> phase_distribution(0.0, 2 * PI);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int l = 0; l < n; ++l) {
                    const int idx = (i * n + j) * n + l;

                    const int kx = (i <= n / 2) ? i : i - n;
                    const int ky = (j <= n / 2) ? j : j - n;
                    const int kz = (l <= n / 2) ? l : l - n;
                    const double k_mag = std::sqrt(static_cast<double>(kx * kx + ky * ky + kz * kz));

                    if (k_mag < 1e-12 || k_mag < k - 1 || k_mag > k) {
                        continue;
                    }

                    const double amp_x = amplitude_distribution(generator);
                    const double amp_y = amplitude_distribution(generator);
                    const double amp_z = amplitude_distribution(generator);
                    const double phase_x = phase_distribution(generator);
                    const double phase_y = phase_distribution(generator);
                    const double phase_z = phase_distribution(generator);

                    const std::complex<double> ax = amp_x * std::complex<double>(std::cos(phase_x), std::sin(phase_x));
                    const std::complex<double> ay = amp_y * std::complex<double>(std::cos(phase_y), std::sin(phase_y));
                    const std::complex<double> az = amp_z * std::complex<double>(std::cos(phase_z), std::sin(phase_z));

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

        fftw_plan plan_x = fftw_plan_dft_3d(n, n, n, fieldx, realx, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_plan plan_y = fftw_plan_dft_3d(n, n, n, fieldy, realy, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_plan plan_z = fftw_plan_dft_3d(n, n, n, fieldz, realz, FFTW_BACKWARD, FFTW_ESTIMATE);

        fftw_execute(plan_x);
        fftw_execute(plan_y);
        fftw_execute(plan_z);

        outx.resize(N);
        outy.resize(N);
        outz.resize(N);
        for (size_t idx = 0; idx < N; ++idx) {
            outx[idx] = realx[idx][0] / static_cast<double>(N);
            outy[idx] = realy[idx][0] / static_cast<double>(N);
            outz[idx] = realz[idx][0] / static_cast<double>(N);
        }

        fftw_destroy_plan(plan_x);
        fftw_destroy_plan(plan_y);
        fftw_destroy_plan(plan_z);
        fftw_free(fieldx);
        fftw_free(fieldy);
        fftw_free(fieldz);
        fftw_free(realx);
        fftw_free(realy);
        fftw_free(realz);
    };

    std::vector<double> solx, soly, solz;
    std::vector<double> compx, compy, compz;

    if (mode == 1) {
        // Single-field Helmholtz mix (Federrath-style)
        build_field(w_sol_mix, w_comp_mix, seed, pertx, perty, pertz);
    } else {
        // Two-field approach: pure sol + pure comp, each normalized to rms=1, then mixed by energy weights
        build_field(1.0, 0.0, seed, solx, soly, solz);
        scale_velocity_to_unity(solx, soly, solz);

        build_field(0.0, 1.0, seed + 1, compx, compy, compz);
        scale_velocity_to_unity(compx, compy, compz);

        pertx.resize(N);
        perty.resize(N);
        pertz.resize(N);
        for (size_t idx = 0; idx < N; ++idx) {
            pertx[idx] = w_sol_mix * solx[idx] + w_comp_mix * compx[idx];
            perty[idx] = w_sol_mix * soly[idx] + w_comp_mix * compy[idx];
            pertz[idx] = w_sol_mix * solz[idx] + w_comp_mix * compz[idx];
        }
    }

}

void scale_velocity_to_unity(std::vector<double>& pertx, std::vector<double>& perty, std::vector<double>& pertz)
{
    const size_t N = pertx.size();
    double current_vrms = 0.0;

    for (size_t i = 0; i < N; ++i) {
        current_vrms += pertx[i] * pertx[i] + perty[i] * perty[i] + pertz[i] * pertz[i];
    }
    current_vrms = std::sqrt(current_vrms / static_cast<double>(N));

    const double scale = (current_vrms > 0.0) ? 1.0 / current_vrms : 1.0;

    for (size_t i = 0; i < N; ++i) {
        pertx[i] *= scale;
        perty[i] *= scale;
        pertz[i] *= scale;
    }
}

void save_to_hdf5(std::vector<double>& pertx, std::vector<double>& perty, std::vector<double>& pertz, int n, int k, unsigned int seed, double f_sol, int mode)
{
    std::stringstream filename;
    filename << "zdrv_n" << n
             << "_k" << k
             << "_seed" << seed
             << "_f" << std::setprecision(6) << f_sol
             << "_mode" << mode
             << ".hdf5";

    const hsize_t dims[3] = {static_cast<hsize_t>(n), static_cast<hsize_t>(n), static_cast<hsize_t>(n)};

    hid_t file = H5Fcreate(filename.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hid_t space = H5Screate_simple(3, dims, nullptr);
    hid_t dset_x = H5Dcreate2(file, "pertx", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t dset_y = H5Dcreate2(file, "perty", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t dset_z = H5Dcreate2(file, "pertz", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    H5Dwrite(dset_x, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pertx.data());
    H5Dwrite(dset_y, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, perty.data());
    H5Dwrite(dset_z, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pertz.data());

    hid_t scalar_space = H5Screate(H5S_SCALAR);
    hid_t attr_n = H5Acreate2(file, "n", H5T_NATIVE_INT, scalar_space, H5P_DEFAULT, H5P_DEFAULT);
    hid_t attr_k = H5Acreate2(file, "k", H5T_NATIVE_INT, scalar_space, H5P_DEFAULT, H5P_DEFAULT);
    hid_t attr_seed = H5Acreate2(file, "seed", H5T_NATIVE_UINT, scalar_space, H5P_DEFAULT, H5P_DEFAULT);
    hid_t attr_fsol = H5Acreate2(file, "f_sol", H5T_NATIVE_DOUBLE, scalar_space, H5P_DEFAULT, H5P_DEFAULT);
    hid_t attr_mode = H5Acreate2(file, "mode", H5T_NATIVE_INT, scalar_space, H5P_DEFAULT, H5P_DEFAULT);

    H5Awrite(attr_n, H5T_NATIVE_INT, &n);
    H5Awrite(attr_k, H5T_NATIVE_INT, &k);
    H5Awrite(attr_seed, H5T_NATIVE_UINT, &seed);
    H5Awrite(attr_fsol, H5T_NATIVE_DOUBLE, &f_sol);
    H5Awrite(attr_mode, H5T_NATIVE_INT, &mode);

    H5Aclose(attr_fsol);
    H5Aclose(attr_seed);
    H5Aclose(attr_k);
    H5Aclose(attr_n);
    H5Aclose(attr_mode);
    H5Sclose(scalar_space);

    H5Dclose(dset_z);
    H5Dclose(dset_y);
    H5Dclose(dset_x);
    H5Sclose(space);
    H5Fclose(file);
}
