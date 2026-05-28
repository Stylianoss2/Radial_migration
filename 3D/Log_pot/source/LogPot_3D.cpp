#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>

#include "LogPot_3D.h"

using std::string;
using std::vector;

// Global axisymmetric potential parameters
double v_c = 0.22;     // kpc/Myr = 220 km/s
double q   = 0.9;      // flattening parameter

// Global spiral perturbation parameters
double m           = 4.0;            //spiral mode
double G           = 4.5e-12;       // kpc^3 / (Msun Myr^2)
double R_o         = 2.5;           // kpc
double R_CR        = 10.0;          // kpc
double Sigma_o     = 1.226e+9;      // Msun / kpc^2
double e_s         = 0.1;           // fractional surface density spiral amplitude
double theta_deg   = 15.0;          // pitch angle in degrees
double theta       = theta_deg * M_PI / 180.0;
double alpha       = m / std::tan(theta);
double omega_p     = 0.022;         // pattern speed in rad/Myr
double pert_strength = 0.1;         // additional numerical scaling factor

// Vertical scale height (Large values make the perturbation almost independent of z)
double spiral_z_scale = 0.3;        // kpc

// Time envelope for the perturbation.
double t_spiral_start     = 2000.0;
double t_spiral_ramp_up   = 500.0;
double t_spiral_hold      = 0.0;
double t_spiral_ramp_down = 500.0;



// Phi = 0.5 * v_c^2 * ln(x^2 + y^2 + z^2/q^2)
double potential_energy(const Vec3& pos)
{
    const double q2 = q * q;
    const double S = pos.x * pos.x + pos.y * pos.y + pos.z * pos.z / q2;

    if (S <= 0.0)
    {
        std::cerr << "Error: logarithmic potential is singular at the origin.\n";
        return 0.0;
    }

    return 0.5 * v_c * v_c * std::log(S);
}

Vec3 LogPot_acc(const Vec3& pos)
{
    const double q2 = q * q;
    const double S = pos.x * pos.x + pos.y * pos.y + pos.z * pos.z / q2;

    if (S <= 0.0)
    {
        std::cerr << "Error: acceleration is singular at the origin.\n";
        return {0.0, 0.0, 0.0};
    }

    const double ax = -v_c * v_c * pos.x / S;
    const double ay = -v_c * v_c * pos.y / S;
    const double az = -v_c * v_c * pos.z / (q2 * S);

    return {ax, ay, az};
}

// ======================================================
// Spiral perturbation
// Phi_sp(R,phi,z,t) = s(t) A(R) Z(z) cos(psi)
// s(t) =  perturbation scaling factor = time profile of simulation
// psi = alpha ln(R/R_CR) + m omega_p t - m phi
// A(R) = perturbation amplitude = pert_strength 2 pi G e_s Sigma(R) / k(R)
// k(R) = alpha / R
// Z(z) = sech^2(z / spiral_z_scale)
// ======================================================

double perturbation_scaling_factor(double sim_time)
{
    const double t_up_end   = t_spiral_start + t_spiral_ramp_up;
    const double t_hold_end = t_up_end + t_spiral_hold;
    const double t_down_end = t_hold_end + t_spiral_ramp_down;

    if (sim_time < t_spiral_start)
        return 0.0;

    if (t_spiral_ramp_up > 0.0 && sim_time < t_up_end)
        return (sim_time - t_spiral_start) / t_spiral_ramp_up;

    if (t_spiral_ramp_up == 0.0 && sim_time < t_hold_end)
        return 1.0;

    if (sim_time < t_hold_end)
        return 1.0;

    if (t_spiral_ramp_down > 0.0 && sim_time < t_down_end)
        return (t_down_end - sim_time) / t_spiral_ramp_down;

    if (t_spiral_ramp_down == 0.0 && sim_time >= t_hold_end)
        return 0.0;

    return 0.0;
}

double spiral_vertical_envelope(double z)
{
    if (spiral_z_scale <= 0.0)
        return 1.0;

    const double u = z / spiral_z_scale;

    // Avoid overflow in cosh for very large z
    if (std::abs(u) > 20.0)
        return 0.0;

    //sech^2(u) = 1 / cosh^2(u) 
    const double c = std::cosh(u);
    return 1.0 / (c * c);
}

double d_spiral_vertical_envelope_dz(double z)
{
    if (spiral_z_scale <= 0.0)
        return 0.0;

    const double Z = spiral_vertical_envelope(z);
    const double u = z / spiral_z_scale;

    //-2sech^2(u)tanh(u)
    return -(2.0 / spiral_z_scale) * Z * std::tanh(u);
}

double spiral_potential_3D(const Vec3& pos, double sim_time)
{
    const double x = pos.x, y = pos.y, z = pos.z;

    const double R2 = x * x + y * y;
    const double R  = std::sqrt(R2);

    if (R < 1e-12)
        return 0.0;

    const double phi = std::atan2(y, x);

    const double Sigma = Sigma_o * std::exp(-R / R_o);
    const double k_abs = std::abs(alpha) / R;
    const double Phi_S = (2.0 * M_PI * G * e_s * Sigma) / k_abs;

    const double Z = spiral_vertical_envelope(z); // Stars away from midplane still affected by perturbations but to a lesser extent
    const double scale = perturbation_scaling_factor(sim_time);

    const double A = Phi_S  * pert_strength * scale;

    const double psi = alpha * std::log(R / R_CR) + m * omega_p * sim_time - m * phi;

    return A * Z * std::cos(psi);
}

Vec3 spiral_acceleration_3D(const Vec3& pos, double sim_time)
{
    const double x = pos.x, y = pos.y, z = pos.z;

    const double R2 = x * x + y * y;
    const double R  = std::sqrt(R2);

    if (R < 1e-12 || pert_strength == 0.0)
        return {0.0, 0.0, 0.0};

    const double phi = std::atan2(y, x);

    const double Sigma = Sigma_o * std::exp(-R / R_o);
    const double k_abs = std::abs(alpha) / R;
    const double Phi_S = (2.0 * M_PI * G * e_s * Sigma) / k_abs;

    const double scale = perturbation_scaling_factor(sim_time);
    if (scale == 0.0)
        return {0.0, 0.0, 0.0};

    const double Z     = spiral_vertical_envelope(z);
    const double dZ_dz = d_spiral_vertical_envelope_dz(z);

    // A(R,t) excludes the vertical envelope.
    const double A = Phi_S * pert_strength * scale;

    // Phi_S ∝ R * exp(-R/R_o) -> dPhi_S/dR = Phi_S * (1/R - 1/R_o).
    const double dPhiS_dR = Phi_S * (1.0 / R - 1.0 / R_o);
    const double dA_dR = pert_strength * scale * dPhiS_dR;

    const double psi = alpha * std::log(R / R_CR) + m * omega_p * sim_time - m * phi;

    const double cpsi = std::cos(psi);
    const double spsi = std::sin(psi);

    const double dpsi_dR = alpha / R;

    // Cylindrical derivatives of Phi_sp.
    const double dPhi_dR   = Z * (dA_dR * cpsi - A * spsi * dpsi_dR);
    const double dPhi_dphi = m * A * Z * spsi;
    const double dPhi_dz   = A * dZ_dz * cpsi;

    // Convert gradient from cylindrical to Cartesian.
    const double dPhi_dx = dPhi_dR * (x / R) + dPhi_dphi * (-y / R2);
    const double dPhi_dy = dPhi_dR * (y / R) + dPhi_dphi * ( x / R2);

    return {-dPhi_dx, -dPhi_dy, -dPhi_dz};
}

double total_potential_energy_3D(const Vec3& pos, double sim_time, bool use_pert)
{
    double Phi = potential_energy(pos);

    if (use_pert)
        Phi += spiral_potential_3D(pos, sim_time);

    return Phi;
}

Vec3 total_acceleration_3D(const Vec3& pos, double sim_time, bool use_pert)
{
    Vec3 a = LogPot_acc(pos);

    if (!use_pert)
        return a;

    Vec3 asp = spiral_acceleration_3D(pos, sim_time);

    return {a.x + asp.x, a.y + asp.y, a.z + asp.z};
}


// Integrators
void Leapfrog_integrator_unperturbed(Vec3& pos, Vec3& vel, const double& dt)
{
    Vec3 r_half
    {
        pos.x + 0.5 * dt * vel.x,
        pos.y + 0.5 * dt * vel.y,
        pos.z + 0.5 * dt * vel.z
    };

    Vec3 a = LogPot_acc(r_half);

    vel.x += dt * a.x;
    vel.y += dt * a.y;
    vel.z += dt * a.z;

    pos.x = r_half.x + 0.5 * dt * vel.x;
    pos.y = r_half.y + 0.5 * dt * vel.y;
    pos.z = r_half.z + 0.5 * dt * vel.z;
}

void Leapfrog_integrator_3D(Vec3& pos, Vec3& vel, const double& dt, double sim_time, bool use_pert)
{
    Vec3 r_half
    {
        pos.x + 0.5 * dt * vel.x,
        pos.y + 0.5 * dt * vel.y,
        pos.z + 0.5 * dt * vel.z
    };

    // Time-centred force for the rotating, explicitly time-dependent perturbation.
    Vec3 a = total_acceleration_3D(r_half, sim_time + 0.5 * dt, use_pert);

    vel.x += dt * a.x;
    vel.y += dt * a.y;
    vel.z += dt * a.z;

    pos.x = r_half.x + 0.5 * dt * vel.x;
    pos.y = r_half.y + 0.5 * dt * vel.y;
    pos.z = r_half.z + 0.5 * dt * vel.z;
}

void leapfrog_step(Vec3& pos, Vec3& vel, double dt)
{
    Leapfrog_integrator_unperturbed(pos, vel, dt);
}


// Diagnostics
double kinetic_energy(const Vec3& vel)
{
    return 0.5 * (vel.x * vel.x + vel.y * vel.y + vel.z * vel.z);
}

double total_energy(const Vec3& pos, const Vec3& vel)
{
    return kinetic_energy(vel) + potential_energy(pos);// No spiral perturbation
}

double total_energy_3D(const Vec3& pos, const Vec3& vel, double sim_time, bool use_pert)
{
    return kinetic_energy(vel) + total_potential_energy_3D(pos, sim_time, use_pert);//with spiral perturbation
}

Vec3 angular_momentum(const Vec3& pos, const Vec3& vel)
{
    return cross_product(pos, vel);
}

double angular_momentum_z(const Vec3& pos, const Vec3& vel)
{
    return pos.x * vel.y - pos.y * vel.x;
}

double jacobi_integral_3D(const Vec3& pos, const Vec3& vel, double sim_time, bool use_pert)
{
    const double E  = total_energy_3D(pos, vel, sim_time, use_pert);
    const double Lz = angular_momentum_z(pos, vel);

    return E - omega_p * Lz;
}


// ======================================================
// Orbit output
// ======================================================

void simulate_trajectory(Vec3& orbit_pos, Vec3& orbit_vel, const std::string& basename, double& total_time, double& dt, bool use_pert)
{
    std::ofstream trajectory_file(basename + "_position.dat");
    std::ofstream parameters_file(basename + "_parameters.dat");
    std::ofstream accel_file     (basename + "_accelerations.dat");
    std::ofstream spiral_file    (basename + "_spiral_profile.dat");

    if (!trajectory_file || !parameters_file || !accel_file || !spiral_file)
    {
        std::cerr << "Error: Cannot open one of the output files.\n";
        return;
    }

    trajectory_file << "# t x y z R phi vx vy vz vR vphi vz_cyl\n";
    parameters_file << "# t Lx Ly Lz KE PE_axis PE_sp PE_total E EJ dE dLz dEJ scale\n";
    accel_file      << "# t ax ay az ax_sp ay_sp az_sp\n";
    spiral_file     << "# t R phi z Phi_sp scale Z\n";

    const int steps = static_cast<int>(total_time / dt);
    const int print_interval = std::max(1, steps / 10);

    const double E0  = total_energy_3D(orbit_pos, orbit_vel, 0.0, use_pert);
    const double Lz0 = angular_momentum_z(orbit_pos, orbit_vel);
    const double EJ0 = jacobi_integral_3D(orbit_pos, orbit_vel, 0.0, use_pert);

    for (int i = 0; i <= steps; i++)
    {
        const double sim_time = static_cast<double>(i) * dt;

        Cyl pos_cyl = cart_to_cyl(orbit_pos);
        Cyl vel_cyl = cart_to_cyl_vel(orbit_pos, orbit_vel);

        Vec3 acc_total = total_acceleration_3D(orbit_pos, sim_time, use_pert);
        Vec3 acc_sp    = use_pert ? spiral_acceleration_3D(orbit_pos, sim_time) : Vec3{0.0, 0.0, 0.0};

        const double KE       = kinetic_energy(orbit_vel);
        const double PE_axis  = potential_energy(orbit_pos);
        const double PE_sp    = use_pert ? spiral_potential_3D(orbit_pos, sim_time) : 0.0;
        const double PE_total = PE_axis + PE_sp;
        const double E        = KE + PE_total;

        Vec3 L = angular_momentum(orbit_pos, orbit_vel);
        const double Lz = L.z;
        const double EJ = E - omega_p * Lz;

        const double dE  = (std::abs(E0)  > 0.0) ? (E  - E0)  / std::abs(E0)  : 0.0;
        const double dLz = (std::abs(Lz0) > 0.0) ? (Lz - Lz0) / std::abs(Lz0) : 0.0;
        const double dEJ = (std::abs(EJ0) > 0.0) ? (EJ - EJ0) / std::abs(EJ0) : 0.0;

        const double scale = use_pert ? perturbation_scaling_factor(sim_time) : 0.0;
        const double Z     = use_pert ? spiral_vertical_envelope(orbit_pos.z)  : 0.0;

        trajectory_file << std::setprecision(10)
                        << sim_time << " "
                        << orbit_pos.x << " " << orbit_pos.y << " " << orbit_pos.z << " "
                        << pos_cyl.R << " " << pos_cyl.phi << " "
                        << orbit_vel.x << " " << orbit_vel.y << " " << orbit_vel.z << " "
                        << vel_cyl.R << " " << vel_cyl.phi << " " << vel_cyl.z << "\n";

        parameters_file << std::setprecision(10)
                        << sim_time << " "
                        << L.x << " " << L.y << " " << L.z << " "
                        << KE << " " << PE_axis << " " << PE_sp << " " << PE_total << " "
                        << E << " " << EJ << " "
                        << dE << " " << dLz << " " << dEJ << " " << scale << "\n";

        accel_file << std::setprecision(10)
                   << sim_time << " "
                   << acc_total.x << " " << acc_total.y << " " << acc_total.z << " "
                   << acc_sp.x << " " << acc_sp.y << " " << acc_sp.z << "\n";

        spiral_file << std::setprecision(10)
                    << sim_time << " "
                    << pos_cyl.R << " " << pos_cyl.phi << " " << orbit_pos.z << " "
                    << PE_sp << " " << scale << " "  << " " << Z << "\n";

        if (i % print_interval == 0)
        {
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "[Orbit: " << basename << "] "
                      << "t = " << sim_time
                      << " | E = " << E
                      << " | Lz = " << Lz
                      << " | EJ = " << EJ
                      << " | dE = " << dE
                      << " | dLz = " << dLz
                      << " | dEJ = " << dEJ
                      << " | scale = " << scale << "\n";
        }

        Leapfrog_integrator_3D(orbit_pos, orbit_vel, dt, sim_time, use_pert);
    }

    trajectory_file.close();
    parameters_file.close();
    accel_file.close();
    spiral_file.close();
}

int main()
{
    double total_time = 4000.0; //Myr
    double dt = 0.01; //Myr
    bool use_pert = true;

    std::cout << "Running 3D logarithmic-potential orbit tests.\n";
    std::cout << "v_c = " << v_c << ", q = " << q << "\n";
    std::cout << "use_pert = " << std::boolalpha << use_pert
              << ", m = " << m
              << ", R_CR = " << R_CR
              << ", omega_p = " << omega_p
              << ", pert_strength = " << pert_strength
              << ", spiral_z_scale = " << spiral_z_scale << "\n";

    vector<Cyl> initial_pos =
    {
        {10.0, 0.0,        0.10},
        {10.0, M_PI / 4.0, 0.20},
        {8.0,  M_PI / 2.0, 0.15},
        {12.0, M_PI,       0.30}
    };

    for (size_t i = 0; i < initial_pos.size(); i++)
    {
        Cyl pos_cyl = initial_pos[i];

        const double v_phi = 0.98 * v_c;
        const double v_R   = 0.02 * v_c;
        const double v_z   = 0.02 * v_c;

        Cyl vel_cyl = {v_R, v_phi, v_z};

        Vec3 pos = cyl_to_cart(pos_cyl);
        Vec3 vel = cyl_to_cart_vel(pos_cyl, vel_cyl);

        std::ostringstream basename;
        basename << "orbit3D_" << i;

        simulate_trajectory(pos, vel, basename.str(), total_time, dt, use_pert);
    }

    std::cout << "Finished all 3D orbit tests.\n";

    return 0;
}
