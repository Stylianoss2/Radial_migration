#pragma once

#include <random>
#include <string>
#include <vector>

#include "VecUtils_3D.h"

// One sampled 3D warm-disc phase-space point.
struct SampleOrbit3D
{
    double R;       // cylindrical radius [kpc]
    double phi;     // azimuth [rad]
    double z;       // height [kpc]

    double v_R;     // radial velocity [kpc/Myr]
    double v_phi;   // azimuthal velocity [kpc/Myr]
    double v_z;     // vertical velocity [kpc/Myr]

    double Lz;      // angular momentum about z-axis [kpc^2/Myr]
    double Rg;      // guiding-centre estimate Lz / v_c [kpc]

    Vec3 pos;       // Cartesian position
    Vec3 vel;       // Cartesian velocity
};

class InitialConditions3D
{
public:
    InitialConditions3D(int N_R_grid = 30000, double R_min = 5.0, double R_max = 15.0);
                        
    // Draw one star from the exponential-disc radial distribution.
    SampleOrbit3D sample_star(std::mt19937& gen);

    // Write N sampled stars to a space-separated .dat file.
    void write_initial_conditions(const std::string& filename, int N_samples, unsigned int seed = 42);
                                  
    // Disc and frequency helpers.
    double Sigma(double R);
    double sigma_R(double R);
    double sigma_z(double R);
    double Omega(double R);
    double kappa(double R);
    double nu(double R);
    double asymmetric_drift(double R);

    // Circular-orbit helpers for the logarithmic potential.
    double E_circ(double R);
    double L_circ(double R);

private:
    void precompute_cdf();
    double sample_R_exponential(std::mt19937& gen);

    std::vector<double> R_vals;
    std::vector<double> CDF_vals;

    int N_grid;
    double R_min;
    double R_max;

    // Disc model parameters.
    double R_o = 2.5;             // exponential surface-density scale length [kpc]
    double R_sigma = 7.5;         // velocity-dispersion scale length [kpc]
    double Sigma_o = 1.226e+9;    // surface-density normalization [Msun/kpc^2]

    // sigma_R0 and sigma_z0 are the dispersions at R_ref.
    double sigma_R0 = 0.02;       // kpc/Myr, about 20 km/s
    double sigma_z0 = 0.006;      // kpc/Myr, about 6 km/s
    double R_ref = 10.0;          // kpc, usually near corotation

    // Potential parameters, matching the logarithmic potential.
    double v_c = 0.22;            // kpc/Myr
    double q = 0.9;               // vertical flattening

    // Rejection limits to avoid extreme initial conditions.
    double max_abs_z = 1.5;       // kpc
    double sigma_cut = 3.0;       // reject |v| > sigma_cut sigma
};
