#include <iostream>
#include <limits>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <sstream>
#include <tuple>
#include <functional>

#include "VecUtils.h"
#include "LogPot.h"

using std::string;
using std::vector;
using std::ofstream;
using GradCandidate = std::tuple<double, double, double>; // [x, y, gradient]
using std::cos;
using std::sin;
using std::tan;
using std::log;

// Global parameters
double m           = 4.0;
double G           = 4.5e-12;       // kpc^3 / (Msun * Myr^2)
double R_o         = 2.5;           // kpc
double R_sigma     = 3 * R_o;       // kpc
double R_CR        = 10;            // kpc
double Sigma_o     = 1.226e+9;      // Msun / kpc^2
double sigma_o     = 0.0;           // σ_R / 2 from Chiba et al. 2021
double e_s         = 0.1;           // old working value 0.3
double theta_deg   = 15;            // old working value 35
double theta       = theta_deg * M_PI / 180.0;
double alpha       = m / std::tan(theta);
double omega_p     = 0.022;         // pattern speed (rad/Myr)
double pert_strength = 0.1;
double v_c         = 0.22;          // kpc/Myr = 220 km/s

static inline Vec2 cyl_vel_to_cart(double vR, double vphi, const Vec2& pos)
{
    double phi = std::atan2(pos.y, pos.x);
    return {
        vR * cos(phi) - vphi * sin(phi),
        vR * sin(phi) + vphi * cos(phi)
    };
}

static inline std::pair<double, double> cart_to_cyl_inertial(const Vec2& pos, const Vec2& vel)
{
    double R = std::hypot(pos.x, pos.y);
    if (R < 1e-12) return {0.0, 0.0};

    double c = pos.x / R;
    double s = pos.y / R;

    double vR   =  c * vel.x + s * vel.y;
    double vphi = -s * vel.x + c * vel.y;

    return {vR, vphi};
}

static inline double guiding_radius_from_L(double Lz)
{
    return Lz / v_c;
}

// Effective potential 
double effective_potential(const Vec2 &cart_pos, const bool &use_pert)
{
    Cyl cyl_pos = cartesian_to_cylindrical(cart_pos);

    // Axisymmetric log potential
    double log_pot_energy = 0.5 * v_c * v_c *
        log(cart_pos.x * cart_pos.x + cart_pos.y * cart_pos.y);

    // Optional spiral perturbation
    double pert_energy = 0.0;
    if (use_pert)
    {
        double kappa = alpha / cyl_pos.R;
        double Sigma = Sigma_o * exp(-cyl_pos.R / R_o);
        double Phi_S = (2.0 * M_PI * G * e_s * Sigma) / kappa;

        // Spiral-phase (no explicit time term in this snapshot version)
        pert_energy = (Phi_S * pert_strength) * cos((alpha * log(cyl_pos.R / R_CR)) - m * cyl_pos.phi);

        // For a bar:
        // pert_energy = (Phi_S * pert_strength) * cos(m * cyl_pos.phi);
    }

    // Rotating-frame centrifugal term
    double cent_term = 0.5 * omega_p * omega_p * (cyl_pos.R * cyl_pos.R);

    return log_pot_energy + pert_energy + cent_term;
}

Vec2 effective_potential_gradient(const Vec2 &cart_pos, const bool &use_pert)
{
    double delta = 1e-8;

    Vec2 dx_pos = {cart_pos.x + delta, cart_pos.y};
    Vec2 dx_neg = {cart_pos.x - delta, cart_pos.y};
    Vec2 dy_pos = {cart_pos.x, cart_pos.y + delta};
    Vec2 dy_neg = {cart_pos.x, cart_pos.y - delta};

    Vec2 grad;
    grad.x = (effective_potential(dx_pos, use_pert) - effective_potential(dx_neg, use_pert)) / (2.0 * delta);
    grad.y = (effective_potential(dy_pos, use_pert) - effective_potential(dy_neg, use_pert)) / (2.0 * delta);

    return grad;
}

// Lagrange point identification
bool is_near_cluster(const vector<GradCandidate>& cluster, const GradCandidate& point, double epsilon)
{
    for (const auto& c : cluster)
    {
        double dx = std::get<0>(c) - std::get<0>(point);
        double dy = std::get<1>(c) - std::get<1>(point);
        if (std::sqrt(dx * dx + dy * dy) < epsilon) return true;
    }
    return false;
}

void scan_Lagrange_points(double xmin, double xmax, double ymin, double ymax, double step, const string& filename, bool use_pert)
{
    std::ofstream fout(filename);
    std::ofstream fout_lagrange("lagrange_points.dat");

    fout << "# x\ty\t|∇Phi_eff|\n";
    fout_lagrange << "# Lagrange points\n";

    vector<GradCandidate> low_grad_points;
    double threshold = 5e-5;
    double epsilon   = 0.5; // clustering radius in kpc

    // Find all gradients below the threshold
    for (double y = ymin; y <= ymax; y += step)
    {
        for (double x = xmin; x <= xmax; x += step)
        {
            Vec2 pos = {x, y};
            Vec2 grad = effective_potential_gradient(pos, use_pert);
            double grad_mag = std::sqrt(grad.x * grad.x + grad.y * grad.y);

            fout << x << " " << y << " " << grad_mag << "\n";

            if (grad_mag < threshold)
                low_grad_points.emplace_back(x, y, grad_mag);
        }
        fout << "\n";
    }

    // Sort by gradient magnitude
    std::sort(low_grad_points.begin(), low_grad_points.end(), [](const GradCandidate& a, const GradCandidate& b)
              { return std::get<2>(a) < std::get<2>(b); });

    // Cluster
    vector<vector<GradCandidate>> clusters;
    for (const auto& pt : low_grad_points)
    {
        bool added = false;
        for (auto& cluster : clusters)
        {
            if (is_near_cluster(cluster, pt, epsilon))
            {
                cluster.push_back(pt);
                added = true;
                break;
            }
        }
        if (!added)
            clusters.push_back({pt});
    }

    // Output one representative (min grad) per cluster
    for (const auto& cluster : clusters)
    {
        auto min_point = *std::min_element(
            cluster.begin(), cluster.end(),
            [](const GradCandidate& a, const GradCandidate& b)
            { return std::get<2>(a) < std::get<2>(b); });

        fout_lagrange << std::get<0>(min_point) << " " << std::get<1>(min_point) << " " << std::get<2>(min_point) << "\n";
    }

    fout.close();
    fout_lagrange.close();
}

// Perturbation scaling profile
double perturbation_scaling_factor(double t, double /*total_time*/)
{
    const double t_start = 2000.0;
    const double t_peak  = 2500.0;
    const double t_end   = 3000.0;

    if (t < t_start || t > t_end) return 0.0;

    if (t < t_peak)  return (t - t_start) / (t_peak - t_start); // ramp up
    return (t_end - t) / (t_end - t_peak);                      // ramp down
}

double perturbation_scaling_factor_with_hold(double t, double /*total_time*/)
{
    const double t_start   = 2000.0;
    const double dur_up    = 200.0;
    const double dur_hold  = 600.0;
    const double dur_down  = 200.0;

    const double t_up_end   = t_start + dur_up;
    const double t_hold_end = t_up_end + dur_hold;
    const double t_down_end = t_hold_end + dur_down;

    if      (t <  t_start)    return 0.0;
    else if (t <  t_up_end)   return (t - t_start) / dur_up;
    else if (t <  t_hold_end) return 1.0;
    else if (t <  t_down_end) return (t_down_end - t) / dur_down;
    else                      return 0.0;
}

// Energies (inertial and rotating-frame effective)
// Not used during ramping version
double total_potential_energy(const Vec2 &cart_pos, const double &sim_time)
{
    Cyl cyl_pos = cartesian_to_cylindrical(cart_pos);

    double log_pot_energy = 0.5 * v_c * v_c *
        log(cart_pos.x * cart_pos.x + cart_pos.y * cart_pos.y);

    double kappa = alpha / cyl_pos.R;
    double Sigma = Sigma_o * exp(-cyl_pos.R / R_o);
    double Phi_S = (2.0 * M_PI * G * e_s * Sigma) / kappa;

    double pert_energy = (Phi_S * pert_strength) * cos(alpha * log(cyl_pos.R / R_CR) + m * omega_p * sim_time - m * cyl_pos.phi);

    // For a bar:
    // double pert_energy = (Phi_S * pert_strength) * cos(m*omega_p*sim_time + m * cyl_pos.phi);

    return log_pot_energy + pert_energy;
}

double potential_energy_unperturbed(const Vec2 &pos)
{
    return 0.5 * v_c * v_c * log(pos.x * pos.x + pos.y * pos.y);
}

// Time-dependent perturbation (ramp)
double total_PE_time(const Vec2 &cart_pos, const double &sim_time, const double &total_time)
{
    if (pert_strength == 0.0)
        return potential_energy_unperturbed(cart_pos);

    Cyl cyl_pos = cartesian_to_cylindrical(cart_pos);

    double log_pot_energy = 0.5 * v_c * v_c * log(cart_pos.x * cart_pos.x + cart_pos.y * cart_pos.y);

    double kappa = alpha / cyl_pos.R;
    double Sigma = Sigma_o * exp(-cyl_pos.R / R_o);
    double Phi_S = (2.0 * M_PI * G * e_s * Sigma) / kappa;

    double scale = perturbation_scaling_factor(sim_time, total_time);

    double pert_energy = (Phi_S * pert_strength * scale) * cos(alpha * log(cyl_pos.R / R_CR) + m * omega_p * sim_time - m * cyl_pos.phi);

    // For a bar:
    // double pert_energy = (Phi_S * pert_strength * scale) * cos(m * omega_p * sim_time - m * cyl_pos.phi);

    return log_pot_energy + pert_energy;
}

// Axisymmetric logarithmic potential acceleration
Vec2 LogPot_acc(const Vec2 &pos)
{
    double R2 = pos.x * pos.x + pos.y * pos.y;
    double dPhi_dx = -(v_c * v_c * pos.x) / R2;
    double dPhi_dy = -(v_c * v_c * pos.y) / R2;
    return {dPhi_dx, dPhi_dy};
}

void Leapfrog_integrator_unperturbed(Vec2 &pos, Vec2 &vel, const double &dt)
{
    // drift
    Vec2 r_half{ pos.x + 0.5 * dt * vel.x, pos.y + 0.5 * dt * vel.y };

    // kick
    Vec2 a = LogPot_acc(r_half);
    vel.x += dt * a.x;
    vel.y += dt * a.y;

    // drift
    pos.x = r_half.x + 0.5 * dt * vel.x;
    pos.y = r_half.y + 0.5 * dt * vel.y;
}

// Finite-difference acceleration from total_PE_time
Vec2 total_acceleration(const Vec2 &cart_pos, const double &sim_time, const double &total_time)
{
    if (pert_strength == 0.0)
    {
        return LogPot_acc(cart_pos);
    }
    double delta = 1e-6;

    Vec2 pos_dx_p{ cart_pos.x + delta, cart_pos.y };
    Vec2 pos_dx_m{ cart_pos.x - delta, cart_pos.y };
    Vec2 pos_dy_p{ cart_pos.x, cart_pos.y + delta };
    Vec2 pos_dy_m{ cart_pos.x, cart_pos.y - delta };

    double dVdx = (total_PE_time(pos_dx_m, sim_time, total_time) - total_PE_time(pos_dx_p, sim_time, total_time)) * 0.5 / delta;
    double dVdy = (total_PE_time(pos_dy_m, sim_time, total_time) - total_PE_time(pos_dy_p, sim_time, total_time)) * 0.5 / delta;

    return { dVdx, dVdy };
}

void Leapfrog_integrator_perturbed(Vec2 &pos, Vec2 &vel, const double &dt, const double &sim_time, const double &total_time)
{
    // drift
    pos.x += 0.5 * dt * vel.x;
    pos.y += 0.5 * dt * vel.y;

    // kick
    Vec2 a = total_acceleration(pos, sim_time, total_time);
    vel.x += dt * a.x;
    vel.y += dt * a.y;

    // drift
    pos.x += 0.5 * dt * vel.x;
    pos.y += 0.5 * dt * vel.y;
}

// Quantity determination
double angular_momentum(const Vec2 &pos, const Vec2 &vel)
{
    return pos.x * vel.y - pos.y * vel.x;
}

double kinetic_energy(const Vec2 &vel)
{
    return 0.5 * (vel.x * vel.x + vel.y * vel.y);
}

// inertial velocity from rotating frame: v_inert = v_rot + Ω_p × r
double kinetic_energy_inertial_frame(const Vec2 &pos, const Vec2 &vel_rotating)
{
    Vec2 v_inert{ vel_rotating.x - omega_p * pos.y, vel_rotating.y + omega_p * pos.x };
    return 0.5 * (v_inert.x * v_inert.x + v_inert.y * v_inert.y);
}

double angular_momentum_inertial_frame(const Vec2 &pos, const Vec2 &vel_rotating)
{
    Vec2 v_inert{ vel_rotating.x - omega_p * pos.y, vel_rotating.y + omega_p * pos.x };
    return pos.x * v_inert.y - pos.y * v_inert.x;
}

double jacobi_integral(const Vec2 &cart_pos, const Vec2 &cart_vel, const double &sim_time, const double &total_time)
{
    if (pert_strength == 0.0)
    {
        double KE = kinetic_energy(cart_vel);
        double PE = potential_energy_unperturbed(cart_pos);
        return KE + PE;
    }

    double KE = kinetic_energy(cart_vel);
    double PE = total_PE_time(cart_pos, sim_time, total_time);
    double Lz = angular_momentum(cart_pos, cart_vel);

    // E_J = E - Ω_p Lz  (rotating-frame constant when Φ is steady in that frame)
    return KE + PE - (omega_p * Lz);
}

double effective_potential_time(const Vec2 &cart_pos, const bool &use_pert, double sim_time, double total_time)
{
    Cyl cyl_pos = cartesian_to_cylindrical(cart_pos);

    double log_pot_energy = 0.5 * v_c * v_c * log(cart_pos.x * cart_pos.x + cart_pos.y * cart_pos.y);

    double pert_energy = 0.0;
    if (use_pert && pert_strength != 0.0)
    {
        double kappa = alpha / cyl_pos.R;
        double Sigma = Sigma_o * exp(-cyl_pos.R / R_o);
        double Phi_S = (2.0 * M_PI * G * e_s * Sigma) / kappa;
        double scale = perturbation_scaling_factor(sim_time, total_time);

        pert_energy = (Phi_S * pert_strength * scale) * cos((alpha * log(cyl_pos.R / R_CR)) + m * omega_p * sim_time - m * cyl_pos.phi);
        // For a bar:
        // pert_energy = (Phi_S * pert_strength * scale) * cos(m * omega_p * sim_time - m * cyl_pos.phi);
    }

    double cent_term = 0.5 * omega_p * omega_p * (cyl_pos.R * cyl_pos.R);

    return log_pot_energy + pert_energy + cent_term;
}

Vec2 effective_potential_gradient_time(const Vec2 &cart_pos, const bool &use_pert, double sim_time, double total_time)
{
    double delta = 1e-8;

    Vec2 dx_pos = {cart_pos.x + delta, cart_pos.y};
    Vec2 dx_neg = {cart_pos.x - delta, cart_pos.y};
    Vec2 dy_pos = {cart_pos.x, cart_pos.y + delta};
    Vec2 dy_neg = {cart_pos.x, cart_pos.y - delta};

    Vec2 grad;
    grad.x = (effective_potential_time(dx_pos, use_pert, sim_time, total_time) -
              effective_potential_time(dx_neg, use_pert, sim_time, total_time)) / (2.0 * delta);

    grad.y = (effective_potential_time(dy_pos, use_pert, sim_time, total_time) - 
              effective_potential_time(dy_neg, use_pert, sim_time, total_time)) / (2.0 * delta);

    return grad;
}

void simulate_trajectory(Vec2 &orbit_pos, Vec2 &orbit_vel, const std::string &basename, double &total_time, double &dt, bool use_pert)
{
    std::ofstream trajectory_file (basename + "_position.dat");
    std::ofstream parameters_file (basename + "_parameters.dat");
    std::ofstream accel_file      (basename + "_accelerations.dat");
    std::ofstream rot_file        (basename + "_position_rot.dat");

    if (!trajectory_file || !parameters_file || !accel_file || !rot_file)
    {
        std::cerr << "Error: Cannot open one of the output files.\n";
        return;
    }

    trajectory_file << "# t\tx\ty\tR\tphi\n";
    parameters_file << "# t\tL_z\tKE\tPE\tE\tE_J\tR_L\n";
    accel_file      << "# t\tax_log\tay_log\ttax_pert\tay_pert\n";
    rot_file        << "# t\tx_rot\ty_rot\n";

    int steps = static_cast<int>(total_time / dt);
    int print_interval = steps / 10;

    for (int i = 0; i <= steps; i++)
    {
        double sim_time = static_cast<double>(i) * dt;
        double scale = perturbation_scaling_factor(sim_time, total_time);
        bool perturb_active = use_pert && pert_strength > 0.0 && scale > 0.0;

        // Accelerations
        Vec2 log_acc = LogPot_acc(orbit_pos);
        Vec2 tot_acc = perturb_active ? total_acceleration(orbit_pos, sim_time, total_time) : log_acc;

        // Energies & Momenta
        double Lz = perturb_active ? angular_momentum_inertial_frame(orbit_pos, orbit_vel) : angular_momentum(orbit_pos, orbit_vel);

        double KE = perturb_active ? kinetic_energy_inertial_frame(orbit_pos, orbit_vel) : kinetic_energy(orbit_vel);

        double PE = perturb_active ? total_PE_time(orbit_pos, sim_time, total_time) : potential_energy_unperturbed(orbit_pos);

        double E  = KE + PE;
        double EJ = perturb_active ? jacobi_integral(orbit_pos, orbit_vel, sim_time, total_time) : E;
        double R_L = Lz / v_c;

        // Rotating frame position
        Cyl pos_cyl = cartesian_to_cylindrical(orbit_pos);
        Cyl pos_rot = pos_cyl;
        pos_rot.phi -= omega_p * sim_time;
        if (pos_rot.phi < 0)            pos_rot.phi += 2.0 * M_PI;
        if (pos_rot.phi >= 2.0 * M_PI)  pos_rot.phi -= 2.0 * M_PI;
        Vec2 pos_rot_cart = cylindrical_to_cartesian(pos_rot);

        // Output
        trajectory_file << sim_time << "\t"<< orbit_pos.x << "\t" << orbit_pos.y << "\t" << pos_cyl.R << "\t" << pos_cyl.phi << "\n";

        parameters_file << sim_time << "\t" << Lz << "\t" << KE << "\t" << PE << "\t" << E << "\t" << EJ << "\t" << R_L << "\n";

        accel_file << sim_time << "\t" << log_acc.x << "\t" << log_acc.y << "\t" << tot_acc.x << "\t" << tot_acc.y << "\n";

        rot_file << sim_time << "\t" << pos_rot_cart.x << "\t" << pos_rot_cart.y << "\n";

        if (i % print_interval == 0)
        {
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "[Orbit: " << basename << "] t = " << sim_time << " | L_z = " << Lz << " | E = " << E << " | E_J = " << EJ << "\n";
        }

        if (perturb_active)
            Leapfrog_integrator_perturbed(orbit_pos, orbit_vel, dt, sim_time, total_time);
        else
            Leapfrog_integrator_unperturbed(orbit_pos, orbit_vel, dt);
    }

    trajectory_file.close();
    parameters_file.close();
    accel_file.close();
    rot_file.close();
}

// Compare force due to logarithmic potential and spiral perturbation
void map_force_ratios(double xmin, double xmax, double ymin, double ymax, double step, double sim_time, const double &total_time,
                      const std::string& filename)
{
    std::ofstream fout(filename);
    fout << "# x\ty\t|F_perturb|/|F_axisym|\n";

    for (double y = ymin; y <= ymax; y += step)
    {
        for (double x = xmin; x <= xmax; x += step)
        {
            Vec2 pos{ x, y };

            Vec2 F_axisym = LogPot_acc(pos);
            Vec2 F_total  = total_acceleration(pos, sim_time, total_time);

            Vec2 F_perturb{ F_total.x - F_axisym.x,  F_total.y - F_axisym.y };

            double mag_axisym = std::hypot(F_axisym.x, F_axisym.y);
            double mag_pert   = std::hypot(F_perturb.x, F_perturb.y);

            double ratio = (mag_axisym > 1e-12) ? (mag_pert / mag_axisym) : 0.0;

            fout << x << "\t" << y << "\t" << ratio << "\n";
        }
        fout << "\n";
    }

    fout.close();
    std::cout << "Wrote force ratio map to: " << filename << "\n";
}

// Map spiral potential to total potential
void map_potential_ratios(double xmin, double xmax, double ymin, double ymax, double step, double sim_time, const double& total_time, const std::string& filename)
{
    std::ofstream fout(filename);
    fout << "# x\ty\t|Phi_pert|/|Phi_tot|\n";

    for (double y = ymin; y <= ymax; y += step)
    {
        for (double x = xmin; x <= xmax; x += step)
        {
            Vec2 pos{ x, y };
            Cyl cyl_pos = cartesian_to_cylindrical(pos);

            double phi_log = 0.5 * v_c * v_c * std::log(pos.x * pos.x + pos.y * pos.y);

            double kappa = alpha / cyl_pos.R;
            double Sigma = Sigma_o * std::exp(-cyl_pos.R / R_o);
            double Phi_S = (2.0 * M_PI * G * e_s * Sigma) / kappa;

            double scale = perturbation_scaling_factor(sim_time, total_time);

            double phi_pert = (Phi_S * pert_strength * scale) * std::cos(alpha * std::log(cyl_pos.R / R_CR) - m * cyl_pos.phi);

            // For a bar:
            // double phi_pert = (Phi_S * pert_strength * scale) * std::cos(m * cyl_pos.phi);

            double phi_tot = phi_log + phi_pert;
            double ratio = (std::abs(phi_tot) > 1e-12) ? std::abs(phi_pert) / std::abs(phi_tot) : 0.0;

            fout << x << '\t' << y << '\t' << ratio << '\n';
        }
        fout << '\n';
    }

    fout.close();
    std::cout << "Wrote potential ratio map to: " << filename << "\n";
}


// JR and random energy
double getJRunperturbed(Vec2 pos0, Vec2 vel, double dt)
{
    double JRback = 0.0;

    // Initial cylindrical coordinates
    double Rold = std::hypot(pos0.x, pos0.y);
    double vold = (pos0.x * vel.x + pos0.y * vel.y) / Rold; // v_R

    // First leapfrog step
    Leapfrog_integrator_unperturbed(pos0, vel, dt);
    double Rmid = std::hypot(pos0.x, pos0.y);
    double vmid = (pos0.x * vel.x + pos0.y * vel.y) / Rmid;

    // Second leapfrog step
    Leapfrog_integrator_unperturbed(pos0, vel, dt);
    double Rnew = std::hypot(pos0.x, pos0.y);
    double vnew = (pos0.x * vel.x + pos0.y * vel.y) / Rnew;

    // Find first extremum (sign change of discrete second derivative)
    while ((Rmid - Rold) * (Rnew - Rmid) > 0.0)
    {
        Rold = Rmid; vold = vmid;
        Rmid = Rnew; vmid = vnew;
        Leapfrog_integrator_unperturbed(pos0, vel, dt);
        Rnew = std::hypot(pos0.x, pos0.y);
        vnew = (pos0.x * vel.x + pos0.y * vel.y) / Rnew;
    }

    // move away from extremum
    Leapfrog_integrator_unperturbed(pos0, vel, dt);
    Rold = Rmid; vold = vmid;
    Rmid = Rnew; vmid = vnew;
    Rnew = std::hypot(pos0.x, pos0.y);
    vnew = (pos0.x * vel.x + pos0.y * vel.y) / Rnew;

    // next extremum
    while ((Rmid - Rold) * (Rnew - Rmid) > 0.0)
    {
        Rold = Rmid; vold = vmid;
        Rmid = Rnew; vmid = vnew;
        Leapfrog_integrator_unperturbed(pos0, vel, dt);
        Rnew = std::hypot(pos0.x, pos0.y);
        vnew = (pos0.x * vel.x + pos0.y * vel.y) / Rnew;

        JRback += (Rnew - Rold) * 0.5 * (vnew + vmid);
    }

    // Continue to next extremum
    Leapfrog_integrator_unperturbed(pos0, vel, dt);
    Rold = Rmid; vold = vmid;
    Rmid = Rnew; vmid = vnew;
    Rnew = std::hypot(pos0.x, pos0.y);
    vnew = (pos0.x * vel.x + pos0.y * vel.y) / Rnew;

    JRback += (Rnew - Rold) * 0.5 * (vnew + vmid);

    while ((Rmid - Rold) * (Rnew - Rmid) > 0.0)
    {
        Rold = Rmid; vold = vmid;
        Rmid = Rnew; vmid = vnew;
        Leapfrog_integrator_unperturbed(pos0, vel, dt);
        Rnew = std::hypot(pos0.x, pos0.y);
        vnew = (pos0.x * vel.x + pos0.y * vel.y) / Rnew;

        JRback += (Rnew - Rold) * 0.5 * (vnew + vmid);
    }

    return JRback / (2.0 * M_PI); // J_R
}

double random_energy(const Vec2 &pos, const Vec2 &vel)
{
    double Lz = angular_momentum(pos, vel);
    double Rg = Lz / v_c;

    double R = std::hypot(pos.x, pos.y);
    double v_R = (pos.x * vel.x + pos.y * vel.y) / R;

    double potential_R  = 0.5 * v_c * v_c * log(R  * R);
    double potential_Rg = 0.5 * v_c * v_c * log(Rg * Rg);

    double eff_pot_R  = potential_R  + (Lz * Lz) / (2.0 * R  * R);
    double eff_pot_Rg = potential_Rg + (Lz * Lz) / (2.0 * Rg * Rg);

    double eR = 0.5 * v_R * v_R + eff_pot_R - eff_pot_Rg;
    return eR;
}

static void integrate_no_io(Vec2 pos, Vec2 vel, double total_time, double dt, bool use_pert, double& Lz0, double& Lz_end, 
                            double& JR0, double& JR_end, double& eR0, double& eR_end, Vec2* pos_end = nullptr, Vec2* vel_end = nullptr)
{
    Lz0 = angular_momentum(pos, vel);
    JR0 = getJRunperturbed(pos, vel, dt);
    eR0 = random_energy(pos, vel);

    int steps = static_cast<int>(total_time / dt);
    for (int i = 0; i < steps; i++)
    {
        const double t = i * dt;
        const bool perturb_active = use_pert && pert_strength > 0.0 && (perturbation_scaling_factor(t, total_time) > 0.0);

        if (perturb_active)
            Leapfrog_integrator_perturbed(pos, vel, dt, t, total_time);
        else
            Leapfrog_integrator_unperturbed(pos, vel, dt);
    }

    Lz_end = angular_momentum(pos, vel);
    JR_end = getJRunperturbed(pos, vel, dt);
    eR_end = random_energy(pos, vel);

    if (pos_end) *pos_end = pos;
    if (vel_end) *vel_end = vel;
}

int main(int argc, char** argv)
{
    bool use_pert = false;
    if (argc > 1)
    {
        std::string arg = argv[1];
        if      (arg == "-pert")     use_pert = true;
        else if (arg == "-no_pert")  use_pert = false;
        else
        {
            std::cerr << "Unknown option: " << arg << "\n";
            return 1;
        }
    }

    // Perturbation and total simulation durations
    double total_time = 5000.0;
    double dt = 0.01;

    std::ofstream fpert("phi_perturbation_time_series.dat");
    std::ofstream frat ("phi_ratio_time_series.dat");
    std::ofstream fstr ("pert_strength_time_series.dat");

    fpert << "# time [Myr]\tPhi_perturb\n";
    frat  << "# time [Myr]\t|Phi_pert| / |Phi_total|\n";
    fstr  << "# time [Myr]\ts(t)=ramp(0..1)\tA(t)=pert_strength*s(t)\n";

    Vec2 pos = cylindrical_to_cartesian({10.0, 0.0});
    for (double t = 0.0; t <= total_time; t += dt)
    {
        Cyl cyl = cartesian_to_cylindrical(pos);

        double kappa = alpha / cyl.R;
        double Sigma = Sigma_o * std::exp(-cyl.R / R_o);
        double Phi_S = (2.0 * M_PI * G * e_s * Sigma) / kappa;
        double scale = perturbation_scaling_factor(t, total_time);
        double phi_pert  = (Phi_S * pert_strength * scale) * std::cos(alpha * std::log(cyl.R / R_CR) + m * omega_p * t - m * cyl.phi);
        double phi_total = 0.5 * v_c * v_c * std::log(cyl.R * cyl.R) + phi_pert;

        fpert << t << "\t" << phi_pert << "\n";
        frat  << t << "\t" << std::fabs(phi_pert) / std::fabs(phi_total) << "\n";
        fstr  << t << "\t" << scale << "\t" << (pert_strength * scale) << "\n";
    }

    fpert.close();
    frat.close();
    fstr.close();

    std::cout << "Wrote phi_perturbation_time_series.dat, " << "phi_ratio_time_series.dat, and " << "pert_strength_time_series.dat\n";

    // --- Phi_eff snapshot ---
    double x_min = -20.0, x_max = 20.0;
    double y_min = -20.0, y_max = 20.0;
    double step  = 0.1;

    std::ofstream fout("phi_eff.dat");
    fout << "# x [kpc]\ty [kpc]\tPhi_eff(x,y)\n";

    for (double y = y_min; y <= y_max; y += step)
    {
        for (double x = x_min; x <= x_max; x += step)
        {
            Vec2 pos_xy{ x, y };
            double snapshot_time = 2500.0; // peak of ramp
            double phi = effective_potential_time(pos_xy, use_pert, snapshot_time, total_time);
            fout << x << "\t" << y << "\t" << phi << "\n";
        }
        fout << "\n";
    }

    fout.close();
    std::cout << "Wrote phi_eff.dat (perturbation " << (use_pert ? "on" : "off") << ")\n";

    // --- Lagrange points ---
    scan_Lagrange_points(x_min, x_max, y_min, y_max, step, "grad_phi_eff.dat", use_pert);
    std::cout << "Wrote grad_phi_eff.dat and lagrange_points.dat\n";

    // --- Batch: DF_initial_conditions -> deltaLz_from_DF.dat ---
    std::ifstream fin("DF_initial_conditions.dat");
    if (!fin)
    {
        std::cerr << "Could not open DF_initial_conditions.dat\n";
        // continue to single-orbit outputs
    }
    else
    {
        std::ofstream out("deltaLz_from_DF.dat");
        out << "# id  R  x  y  vR  vphi  Lz0  Lz_end  dLz  JR0  JR_end  dJR  eR0  eR_end  deR  vR1\n";

        std::string line;
        int id = 0;
        std::size_t N_rows = 0;
        double sum_dLz = 0.0, sum_dJR = 0.0;

        while (std::getline(fin, line))
        {
            if (line.empty() || line[0] == '#') continue;

            double R, L, vR, vphi, g1, sigma_corr, g2, x, y;
            std::istringstream iss(line);
            if (!(iss >> R >> L >> vR >> vphi >> g1 >> sigma_corr >> g2 >> x >> y))
                continue;

            Vec2 pos0{ x, y };
            Vec2 vel0 = cyl_vel_to_cart(vR, vphi, pos0);

            Vec2 pos_end, vel_end;
            double Lz0 = 0.0, Lz_end = 0.0;
            double JR0 = 0.0, JR_end = 0.0;
            double eR0 = 0.0, eR_end = 0.0;

            integrate_no_io(pos0, vel0, total_time, dt, use_pert, Lz0, Lz_end, JR0, JR_end, eR0, eR_end, &pos_end, &vel_end);

            double dLz = (Lz_end - Lz0);
            double dJR = (JR_end - JR0);
            double deR = (eR_end - eR0);

            double vR1 = cart_to_cyl_inertial(pos_end, vel_end).first;

            out << id << ' ' << R << ' ' << x << ' ' << y << ' ' << vR << ' ' << vphi << ' ' << Lz0 << ' ' << Lz_end << ' ' << dLz << ' '
                << JR0 << ' ' << JR_end << ' ' << dJR << ' ' << eR0 << ' ' << eR_end << ' ' << deR << ' ' << vR1 << '\n';

            sum_dLz += dLz;
            sum_dJR += dJR;
            ++id;
            ++N_rows;
        }

        out.close();
        fin.close();

        double mean_dLz = (N_rows > 0) ? (sum_dLz / N_rows) : 0.0;
        std::cout << "Average ΔLz = " << mean_dLz << " (from " << N_rows << " stars)\n";
        std::cout << "Wrote deltaLz_from_DF.dat (combined Lz, JR, eR columns)\n";
    }

    // Collect samples from this run for RMS in a fixed CR window
    std::ifstream fin2("deltaLz_from_DF.dat");
    if (!fin2)
    {
        std::cerr << "Could not open deltaLz_from_DF.dat for RMS calc\n";
    }
    else
    {
        double Lz0, dLz;
        double sum2 = 0.0;
        int N = 0;

        std::string line;
        while (std::getline(fin2, line))
        {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream iss(line);

            int id;
            double R, x, y, vR, vphi, Lz_end, JR0, JR_end, dJR, eR0, eR_end, deR, vR1;
            if (!(iss >> id >> R >> x >> y >> vR >> vphi >> Lz0 >> Lz_end >> dLz >> JR0 >> JR_end >> dJR >> eR0 >> eR_end >> deR >> vR1))
                continue;

            if (Lz0 >= 2.0 && Lz0 <= 2.4) // fixed CR window
            {
                sum2 += dLz * dLz;
                ++N;
            }
        }

        if (N > 0)
        {
            double rms = std::sqrt(sum2 / double(N));
            std::cout << "RMS(ΔLz) in CR window [2.0, 2.4] = " << rms << " (N=" << N << ")\n";
        }
        else
        {
            std::cout << "No stars in CR window [2.0, 2.4]\n";
        }
    }

    // Initialize a few orbits close to corotation
    std::vector<Vec2> initial_pos =
    {
        { 10.2,  0.0 },
        { -8.8, -4.5 },
        { -4.6,  8.8 },
        {  8.6,  4.4 }
    };

    for (size_t i = 0; i < initial_pos.size(); i++)
    {
        Vec2 p = initial_pos[i];
        double phi = std::atan2(p.y, p.x);

        double v_phi = 0.98 * v_c;
        double v_R   = 0.02 * v_c;

        Vec2 v;
        v.x = v_R * std::cos(phi) - v_phi * std::sin(phi);
        v.y = v_R * std::sin(phi) + v_phi * std::cos(phi);

        std::ostringstream basename;
        basename << "orbit_L" << i;
        simulate_trajectory(p, v, basename.str(), total_time, dt, use_pert);
    }

    return 0;
}
