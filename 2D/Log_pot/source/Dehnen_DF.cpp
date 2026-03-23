#include "Dehnen_DF.h"
#include "VecUtils.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <stdexcept>
#include <cmath>
#include <algorithm>


double WarmDiskDF::Sigma(double R) const 
{
    return Sigma_o * std::exp(-R / R_o);
}

void WarmDiskDF::precompute_cdf() 
{
    if (N_grid < 2) 
    {
        throw std::runtime_error("WarmDiskDF: N_grid must be >= 2");
    }
    if (!(R_max > R_min)) 
    {
        throw std::runtime_error("WarmDiskDF: require R_max > R_min");
    }
    R_vals.resize(N_grid);
    CDF_vals.resize(N_grid);
    double dR = (R_max - R_min) / (N_grid - 1);
    double total = 0.0;

    for (int i = 0; i < N_grid; i++) 
    {
        double R = R_min + i * dR;
        R_vals[i] = R;

        // Shu: weight ∝ R * Σ(R)
        double weight = R * Sigma(R);
        total += weight;
        CDF_vals[i] = total;
    }
    for (double& v : CDF_vals) v /= total;
}

double WarmDiskDF::sample_R(std::mt19937& gen) 
{
    if (CDF_vals.empty() || R_vals.size() != CDF_vals.size()) 
    {
        throw std::runtime_error("Distribution values not initialized correclty");
    }
    std::uniform_real_distribution<> dist(0.0, 1.0);
    double u = dist(gen);

    // Binary search for u in CDF
    auto it = std::lower_bound(CDF_vals.begin(), CDF_vals.end(), u);
    int idx = static_cast<int>(std::distance(CDF_vals.begin(), it));
    if (idx < 1) idx = 1;
    if (idx >= N_grid) idx = N_grid - 1;

    double R1 = R_vals[idx - 1];
    double R2 = R_vals[idx];
    double C1 = CDF_vals[idx - 1];
    double C2 = CDF_vals[idx];
    const double denom = (C2 - C1);
    if (denom <= 0.0) 
    {
        // fallback: return left edge of the interval
        return R1;
    }
    double t = (u - C1) / denom;
    return R1 + t * (R2 - R1);
}

double WarmDiskDF::E_circ(double R) const 
{
    return 0.5 * v_c * v_c + 0.5 * v_c * v_c * std::log(R * R);
}

double WarmDiskDF::L_circ(double R) const 
{
    return R * v_c;
}

double WarmDiskDF::Omega(double R) const 
{
    return v_c / R;
}

// Constructor to initialize values
WarmDiskDF::WarmDiskDF(int N_grid_, double R_min_, double R_max_) : N_grid(N_grid_), R_min(R_min_), R_max(R_max_)
{
    precompute_cdf();
}

double WarmDiskDF::sigma_R(double R) const
{
    return sigma_o * std::exp(-R / R_sigma);
}

SampleOrbit WarmDiskDF::sample_orbit(std::mt19937& gen)
{
    /*local epicyclic sampler 
    v_R ~ N(0, σ_R)
    v_φ ~ N(<v_φ>, σ_φ) with  σ_φ^2 ≈ (κ^2 / 4Ω^2) σ_R^2
    <v_φ> = v_c - v_a,  v_a ≈ (σ_R^2 / 2 v_c) * [ -d ln Σ/d ln R - d ln σ_R^2/d ln R - (1 - σ_φ^2/σ_R^2) ]   */

    //pick radius from your Σ-weighted CDF
    double R = sample_R(gen);
    if (!std::isfinite(R) || R <= 0.0) R = std::max(1e-6, R);


    const double vc    = v_c;         
    const double Om    = Omega(R);    
    const double kap   = kappa(R);   

    const double sR    = sigma_R(R);                      // target σ_R(R)
    const double sphi  = sR * (kap / (2.0 * Om));         // epicyclic: σ_φ = σ_R * κ / (2Ω)
    const double sR2   = sR * sR;
    const double sphi2 = sphi * sphi;

    //gradients for asymmetric drift 
    // Σ(R) = Σ0 e^{-R/R_o}          -> d ln Σ / d ln R = - R / R_o
    // σ_R(R) = σ0 e^{-R/R_σ}       -> d ln σ_R^2 / d ln R = - 2 R / R_σ
    const double dlnSigma_dlnR   = - R / R_o;
    const double dlnsR2_dlnR     = - 2.0 * R / R_sigma;

    //asymmetric drift 
    const double bracket = -dlnSigma_dlnR - dlnsR2_dlnR - (1.0 - sphi2 / sR2);
    double v_a = 0.5 * sR2 / std::max(1e-12, vc) * bracket;

    const double v_a_max = 0.5 * vc;                 
    if (!std::isfinite(v_a)) v_a = 0.0;
    v_a = std::clamp(v_a, 0.0, v_a_max);

    const double mean_vphi = vc - v_a;

    std::normal_distribution<> nR(0.0,  sR);
    std::normal_distribution<> nP(mean_vphi, sphi);

    double v_R   = nR(gen);
    double v_phi = nP(gen);

    double L = R * v_phi;

    return { R, L, v_R, v_phi };
}

double WarmDiskDF::kappa(double R) const
{
    return std::sqrt(2.0) * Omega(R);
}

/*
int main()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    WarmDiskDF df;

    const int N_samples = 30000;

    // Open .dat file (space-separated)
    std::ofstream DF("DF_initial_conditions.dat");
    if (!DF) 
    {
        std::cerr << "Error: could not open DF_initial_conditions.dat for writing.\n";
        return 1;
    }

    DF << std::setprecision(6);
    DF << "# R  L  v_R  v_phi   x  y\n";
    std::uniform_real_distribution<double> uniform_phi(0.0, 2.0 * M_PI); //generate random φ values for a specific R coordinate

    for (int i = 0; i < N_samples; i++)
    {
        SampleOrbit orb = df.sample_orbit(gen);
        double phi = uniform_phi(gen);
        Cyl cyl_coords{orb.R, phi};
        Vec2 pos = cylindrical_to_cartesian(cyl_coords);

        DF << orb.R << " " << orb.L << " " << orb.v_R << " " << orb.v_phi <<  " "  << pos.x << " " << pos.y << "\n"; 
    }

    DF.close();
    std::cout << "Generated " << N_samples << " orbits based on Shu Distribution Function\n";
    return 0;
}
  */  


