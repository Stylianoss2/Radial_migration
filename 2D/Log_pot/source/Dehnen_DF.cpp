#include "Dehnen_DF.h"
#include "LogPot.h"
#include "VecUtils.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <numeric>
#include <vector>
#include <cmath>
#include <algorithm>

double WarmDiskDF::Sigma(double R) const 
{
    return Sigma_o * std::exp(-R / R_o);
}


//function for Σ prime
double WarmDiskDF::Sigma_prime(double R) const 
{
    // Prevent R from going beyond grid bounds
    if (R <= R_vals.front()) 
    {
        return Sigma_prime_vals.front();
    }
    if (R >= R_vals.back()) 
    {
        return Sigma_prime_vals.back();
    }

    auto iterator = std::lower_bound(R_vals.begin(), R_vals.end(), R); //find the first element of R in the grid
    int idx = std::distance(R_vals.begin(), iterator); 
    
    double R1 = R_vals[idx - 1];
    double R2 = R_vals[idx];
    double S1 = Sigma_prime_vals[idx - 1];
    double S2 = Sigma_prime_vals[idx];
    double t = (R - R1) / (R2 - R1);
    return S1 + t * (S2 - S1);
}

//function for σ_R prime
double WarmDiskDF::sigmaR_prime(double R) const
{
    if (R <= R_vals.front()) 
    {
        return sigmaR_prime_vals.front();
    }
    if (R >= R_vals.back()) 
    {
        return sigmaR_prime_vals.back();
    }

    auto iterator = std::lower_bound(R_vals.begin(), R_vals.end(), R);
    int idx = std::distance(R_vals.begin(), iterator);

    double R1 = R_vals[idx - 1];
    double R2 = R_vals[idx];
    double S1 = sigmaR_prime_vals[idx - 1];
    double S2 = sigmaR_prime_vals[idx];
    double t = (R - R1) / (R2 - R1);
    return S1 + t * (S2 - S1);
}

// Step 1 in Dehnen-> Find R using a probability distribution
void WarmDiskDF::precompute_cdf() 
{
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
    double t = (u - C1) / (C2 - C1);
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

void WarmDiskDF::update_parameters(const std::vector<double>& Sigma_f, const std::vector<double>& sigmaR_f) 
{
    for (int i = 0; i < N_grid; i++) 
    {
        if (Sigma_f[i] > 0.0) 
        {
            Sigma_prime_vals[i] *= Sigma(R_vals[i]) / Sigma_f[i];
        }
        if (sigmaR_f[i] > 0.0) 
        {
            sigmaR_prime_vals[i] *= sigma_R(R_vals[i]) / sigmaR_f[i];
        }
    }
    precompute_cdf();  
}

// Constructor
WarmDiskDF::WarmDiskDF(int N_grid_, double R_min_, double R_max_) : N_grid(N_grid_), R_min(R_min_), R_max(R_max_)
{
    Sigma_prime_vals.resize(N_grid);
    sigmaR_prime_vals.resize(N_grid);

    double dR = (R_max - R_min) / (N_grid - 1);
    for (int i = 0; i < N_grid; ++i)
    {
        double R = R_min + i * dR;
        Sigma_prime_vals[i] = Sigma(R);
        sigmaR_prime_vals[i] = sigma_R(R);
    }

    precompute_cdf();
}

double WarmDiskDF::sigma_R(double R) const
{
    return sigma_o * std::exp(-R / R_sigma);
}

SampleOrbit WarmDiskDF::sample_orbit(std::mt19937& gen)
{
    // local epicyclic sampler (Method A)
    // v_R ~ N(0, σ_R)
    // v_φ ~ N(<v_φ>, σ_φ) with  σ_φ^2 ≈ (κ^2 / 4Ω^2) σ_R^2
    // and   <v_φ> = v_c - v_a,  v_a ≈ (σ_R^2 / 2 v_c) * [ -d ln Σ/d ln R - d ln σ_R^2/d ln R - (1 - σ_φ^2/σ_R^2) ]

    // 1) pick radius from your Σ-weighted CDF
    double R = sample_R(gen);
    if (!std::isfinite(R) || R <= 0.0) R = std::max(1e-6, R);

    // 2) local circular quantities
    const double vc    = v_c;         // constant in your log potential
    const double Om    = Omega(R);    // = v_c / R
    const double kap   = kappa(R);    // = sqrt(2) * Omega(R) in your code

    // 3) dispersions
    const double sR    = sigma_R(R);                      // target σ_R(R)
    const double sphi  = sR * (kap / (2.0 * Om));         // epicyclic: σ_φ = σ_R * κ / (2Ω)
    const double sR2   = sR * sR;
    const double sphi2 = sphi * sphi;

    // 4) gradients for asymmetric drift (your analytic profiles)
    // Σ(R) = Σ0 e^{-R/R_o}          -> d ln Σ / d ln R = - R / R_o
    // σ_R(R) = σ0 e^{-R/R_σ}       -> d ln σ_R^2 / d ln R = - 2 R / R_σ
    const double dlnSigma_dlnR   = - R / R_o;
    const double dlnsR2_dlnR     = - 2.0 * R / R_sigma;

    // 5) asymmetric drift (Jeans/epicycle form; sign chosen so v_a>0 for falling Σ,σ_R)
    const double bracket = -dlnSigma_dlnR - dlnsR2_dlnR - (1.0 - sphi2 / sR2);
    double v_a = 0.5 * sR2 / std::max(1e-12, vc) * bracket;

    // keep things sane if parameters are extreme
    const double v_a_max = 0.5 * vc;                 // don’t let the drift get crazy
    if (!std::isfinite(v_a)) v_a = 0.0;
    v_a = std::clamp(v_a, 0.0, v_a_max);

    const double mean_vphi = vc - v_a;

    // 6) draw velocities
    std::normal_distribution<> nR(0.0,  sR);
    std::normal_distribution<> nP(mean_vphi, sphi);

    double v_R   = nR(gen);
    double v_phi = nP(gen);

    // 7) angular momentum
    double L = R * v_phi;

    // Dehnen book-keeping factors (not used here)
    double g1 = 1.0, sigma_corr = 0.0, g2 = 1.0;

    return { R, L, v_R, v_phi, g1, sigma_corr, g2 };
}

double WarmDiskDF::kappa(double R) const
{
    return std::sqrt(2.0) * Omega(R);
}

struct OrbitState { Vec2 pos; Vec2 vel; };  
struct PhaseSample { Vec2 pos; Vec2 vel; double R, phi, vR, vphi; };

struct RTable 
{
    // One full radial period sampled during the orbit integration
    std::vector<double> t;     // monotonically increasing, t[0]=0, t.back()≈TR
    std::vector<double> R;     // R(t)
    std::vector<double> Rdot;  // dR/dt (if you didn't store this, you can finite-difference R)
};

static inline void interp_R_Rdot_at(const RTable& tab, double tquery, double& R_out, double& Rdot_out)
{
    // wrap t into [0, t.back())
    double TR = tab.t.back();
    if (TR <= 0) { R_out = tab.R.front(); Rdot_out = tab.Rdot.front(); return; }
    double u = fmod(tquery, TR);
    if (u < 0) u += TR;

    // binary search for interval [i, i+1] with t[i] <= u < t[i+1]
    auto it = std::upper_bound(tab.t.begin(), tab.t.end(), u);
    size_t i = (it == tab.t.begin()) ? 0 : size_t(it - tab.t.begin() - 1);
    size_t j = std::min(i + 1, tab.t.size() - 1);

    double t0 = tab.t[i], t1 = tab.t[j];
    double w  = (t1 > t0) ? (u - t0) / (t1 - t0) : 0.0;

    R_out    = (1.0 - w) * tab.R[i]    + w * tab.R[j];
    Rdot_out = (1.0 - w) * tab.Rdot[i] + w * tab.Rdot[j];
}

static inline int draw_nsamples(double target_mean, std::mt19937_64& rng)
{
    double nf = std::floor(target_mean);
    double p  = target_mean - nf;  // probability of using ceil
    std::uniform_real_distribution<double> U01(0.0, 1.0);
    return (U01(rng) < p) ? int(nf) + 1 : int(nf);
}

std::vector<PhaseSample>
sample_phase_points(const RTable& tab, double TR, double L, int Norb, double g1, double g2, uint64_t seed)
{
    std::mt19937_64 rng(seed);
    const double Ntarget = g1 * g2 * Norb;
    const int Nsam = draw_nsamples(Ntarget, rng);

    std::uniform_real_distribution<double> Uphi(0.0, 2.0 * M_PI);
    std::uniform_real_distribution<double> Ut(0.0, TR);

    std::vector<PhaseSample> out;
    out.reserve(Nsam);

    for (int k = 0; k < Nsam; k++) 
    {
        double phi = Uphi(rng);
        double t   = Ut(rng);

        double R, vR;
        interp_R_Rdot_at(tab, t, R, vR);

        // Guard against pathological R
        if (R <= 0) R = 1e-12;

        // phi-dot and tangential speed
        const double phidot = L / (R * R);
        const double vphi   = L / R; // same thing, just more numerically stable

        // Convert to Cartesian
        const double c = std::cos(phi), s = std::sin(phi);
        Vec2 pos{ R * c, R * s };
        Vec2 vel{
            // v_R * e_R + v_phi * e_phi
            vR * c - vphi * s,
            vR * s + vphi * c
        };

        out.push_back(PhaseSample{pos, vel, R, phi, vR, vphi});
    }

    return out;
}

double compute_radial_period(const OrbitState& initial, double dt, double max_time = 10000.0)
{
    Vec2 pos = initial.pos;
    Vec2 vel = initial.vel;

    double last_peri_time = std::numeric_limits<double>::quiet_NaN();
    double time = 0.0;

    // keep a tiny window for quadratic interpolation
    double R_prev2 = std::numeric_limits<double>::quiet_NaN();
    double R_prev1 = std::numeric_limits<double>::quiet_NaN();
    double R_curr  = std::numeric_limits<double>::quiet_NaN();

    const int max_steps = static_cast<int>(max_time / dt);
    int since_last_ext = 1e+9;             // big number to start
    const int min_sep_samples = 10;       // guard against noise; tune as needed

    for (int i = 0; i < max_steps; i++) 
    {
        // radius now
        double R = std::sqrt(pos.x*pos.x + pos.y*pos.y);

        // shift the window
        R_prev2 = R_prev1;
        R_prev1 = R_curr;
        R_curr  = R;

        // integrate one step
        Leapfrog_integrator_unperturbed(pos, vel, dt);
        time += dt;
        since_last_ext++;

        // need three filled values before testing
        if (!std::isnan(R_prev2) && !std::isnan(R_prev1)) 
        {
            // local minimum at middle sample?
            if (R_prev1 < R_prev2 && R_prev1 < R_curr && since_last_ext >= min_sep_samples) 
            {
                // Quadratic interpolation around (t-2dt, t-dt, t)
                // Fit a parabola through ( -2, R_prev2 ), ( -1, R_prev1 ), ( 0, R_curr ) in units of dt
                double y0 = R_prev2, y1 = R_prev1, y2 = R_curr;
                // time offset of vertex relative to the middle point (-1*dt)
                double denom = (y0 - 2*y1 + y2);
                double tau = 0.0; // in units of dt
                if (std::abs(denom) > 0.0) 
                {
                    // standard 3-point quadratic extremum: t* ≈ (y0 - y2)/(2(y0 - 2y1 + y2))
                    tau = (y0 - y2) / (2.0 * denom);
                }
                // vertex time relative to current 'time' is at t* = time - dt + tau*dt
                double t_peri = time - dt + tau*dt;

                if (std::isnan(last_peri_time)) 
                {
                    last_peri_time = t_peri;          // first pericentre: store and continue
                    since_last_ext = 0;
                } else {
                    double period = t_peri - last_peri_time;
                    if (period > 0.0) return period;  // second pericentre: we're done
                    // if negative (shouldn't happen), ignore and keep going
                    last_peri_time = t_peri;
                    since_last_ext = 0;
                }
            }
        }
    }

    return -1.0; // no period found
}

double compute_g2_from_sample(double R_sample, double v_R, double v_phi, std::function<double(double)> kappa_func) 
{
    // Construct orbit state in x-y plane
    OrbitState orbit;
    orbit.pos = {R_sample, 0.0};  // R along x
    orbit.vel = {v_R, v_phi};     // v_R in x, v_phi in y

    double dt = 0.1; // Adjust as needed
    double TR = compute_radial_period(orbit, dt);

    if (TR <= 0.0)
        return 1.0; // fallback

    double omega_R = 2.0 * M_PI / TR;
    double kappa_val = kappa_func(R_sample);

    return kappa_val / omega_R;
}
/*
int main()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    WarmDiskDF df;

    const int N_samples = 30000;
    const bool compute_g2_flag = true; 

    // Open .dat file (space-separated)
    std::ofstream DF("DF_initial_conditions.dat");
    if (!DF) 
    {
        std::cerr << "Error: could not open DF_initial_conditions.dat for writing.\n";
        return 1;
    }

    DF << std::setprecision(6);
    DF << "# R  L  v_R  v_phi  g1  sigma_corr  g2  x  y\n";
    std::uniform_real_distribution<double> uniform_phi(0.0, 2.0 * M_PI); //generate random φ values for a specific R coordinate

    for (int i = 0; i < N_samples; i++)
    {
        SampleOrbit orb = df.sample_orbit(gen);
        double phi = uniform_phi(gen);
        Cyl cyl_coords{orb.R, phi};
        Vec2 pos = cylindrical_to_cartesian(cyl_coords);

        double g2 = 1.0;
        if (compute_g2_flag) 
        {
            g2 = compute_g2_from_sample( orb.R, orb.v_R, orb.v_phi, [&](double R){ return df.kappa(R); } );
        }

            DF << orb.R << " " << orb.L << " " << orb.v_R << " " << orb.v_phi << " " << orb.g1 << " " << orb.sigma_corr << " "
               << g2 << " " << pos.x << " " << pos.y << "\n"; 
    }

    DF.close();
    std::cout << "Generated " << N_samples << " orbits based on Dehnen Distribution Function\n";
    return 0;
}
*/

