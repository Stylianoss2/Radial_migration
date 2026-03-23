#pragma once
#include <random>
#include <vector>


struct SampleOrbit 
{
    double R; double L; double v_R; double v_phi;
};

class WarmDiskDF 
{
public:
    WarmDiskDF(int N_R_grid = 30000, double R_min = 2.0, double R_max = 20.0);  // N_R_grid is the resolution for CDF. It is not equal to the number of samples
    double sample_R(std::mt19937& gen);
    SampleOrbit sample_orbit(std::mt19937& gen);
    
    double Sigma(double R) const;
    double sigma_R(double R) const;
    double Omega(double R) const;
    double E_circ(double R) const;
    double L_circ(double R) const;
    double kappa(double R) const;

private:
    void precompute_cdf();

    std::vector<double> R_vals;
    std::vector<double> CDF_vals;

    int N_grid;
    double R_min, R_max;

    // Disc model parameters
    double R_o = 2.5;
    double R_sigma = 3.0 * R_o;
    double Sigma_o = 1.226e+9;
    double sigma_o = 0.02;
    double v_c = 0.22;
};