#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <string>

#include "VecUtils_3D.h"

using std::string;
using std::vector;

// Global parameters
double v_c = 0.22;     // kpc/Myr = 220 km/s
double q   = 0.9;      // flattening parameter

// ======================================================
// 3D flattened logarithmic potential
//
// Phi = 0.5 v_c^2 ln(x^2 + y^2 + z^2/q^2)
// ======================================================

double potential_energy(const Vec3& pos)
{
    double q2 = q * q;

    double log_term = pos.x * pos.x + pos.y * pos.y + pos.z * pos.z / q2;

    if (log_term <= 0.0)
    {
        std::cerr << "Error: logarithmic potential is singular at the origin.\n";
        return 0.0;
    }

    return 0.5 * v_c * v_c * std::log(log_term);
}

Vec3 LogPot_acc(const Vec3& pos)
{
    double q2 = q * q;

    double S = pos.x * pos.x + pos.y * pos.y + pos.z * pos.z / q2;

    if (S <= 0.0)
    {
        std::cerr << "Error: acceleration is singular at the origin.\n";
        return {0.0, 0.0, 0.0};
    }

    double ax = -v_c * v_c * pos.x / S;
    double ay = -v_c * v_c * pos.y / S;
    double az = -v_c * v_c * pos.z / (q2 * S);

    return {ax, ay, az};
}

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

// Diagnostics
double kinetic_energy(const Vec3& vel)
{
    return 0.5 * (vel.x * vel.x + vel.y * vel.y + vel.z * vel.z);
}

double total_energy(const Vec3& pos, const Vec3& vel)
{
    return kinetic_energy(vel) + potential_energy(pos);
}

Vec3 angular_momentum(const Vec3& pos, const Vec3& vel)
{
    return cross_product(pos, vel);
}

double angular_momentum_z(const Vec3& pos, const Vec3& vel)
{
    return pos.x * vel.y - pos.y * vel.x;
}

double circular_speed(double R)
{
    // For Phi = 0.5 v_c^2 ln(R^2 + z^2/q^2), in the midplane z = 0, circular speed is constant.
    return v_c;
}

// ======================================================
// Orbit output
// ======================================================

void simulate_trajectory(Vec3& orbit_pos, Vec3& orbit_vel, const std::string& basename, double& total_time, double& dt)
{
    std::ofstream trajectory_file(basename + "_position.dat");
    std::ofstream parameters_file(basename + "_parameters.dat");
    std::ofstream accel_file     (basename + "_accelerations.dat");

    if (!trajectory_file || !parameters_file || !accel_file)
    {
        std::cerr << "Error: Cannot open one of the output files.\n";
        return;
    }

    trajectory_file << "# t x y z R phi vx vy vz vR vphi vz_cyl\n";
    parameters_file << "# t Lx Ly Lz KE PE E dE dLz\n";
    accel_file      << "# t ax ay az\n";

    int steps = static_cast<int>(total_time / dt);
    int print_interval = std::max(1, steps / 10);

    double E0  = total_energy(orbit_pos, orbit_vel);
    double Lz0 = angular_momentum_z(orbit_pos, orbit_vel);

    for (int i = 0; i <= steps; i++)
    {
        double sim_time = static_cast<double>(i) * dt;

        Cyl pos_cyl = cart_to_cyl(orbit_pos);
        Cyl vel_cyl = cart_to_cyl_vel(orbit_pos, orbit_vel);

        Vec3 acc = LogPot_acc(orbit_pos);

        double KE = kinetic_energy(orbit_vel);
        double PE = potential_energy(orbit_pos);
        double E  = KE + PE;

        Vec3 L = angular_momentum(orbit_pos, orbit_vel);
        double Lz = L.z;

        double dE = (E - E0) / std::abs(E0);
        double dLz = (Lz - Lz0) / std::abs(Lz0);

        trajectory_file << std::setprecision(10) << sim_time << " " << orbit_pos.x << " " << orbit_pos.y << " " << orbit_pos.z << " "
                        << pos_cyl.R << " " << pos_cyl.phi << " " << orbit_vel.x << " " << orbit_vel.y << " " << orbit_vel.z << " "
                        << vel_cyl.R << " " << vel_cyl.phi << " " << vel_cyl.z << "\n";

        parameters_file << std::setprecision(10) << sim_time << " " << L.x << " " << L.y << " " << L.z << " " << KE << " " << PE << " " << E << " "
                        << dE << " " << dLz << "\n";

        accel_file << std::setprecision(10) << sim_time << " " << acc.x << " " << acc.y << " " << acc.z << "\n";

        if (i % print_interval == 0)
        {
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "[Orbit: " << basename << "] " << "t = " << sim_time << " | E = " << E << " | Lz = " << Lz << " | dE = " << dE << " | dLz = " << dLz << "\n";            
        }
        Leapfrog_integrator_unperturbed(orbit_pos, orbit_vel, dt);
    }

    trajectory_file.close();
    parameters_file.close();
    accel_file.close();
}

int main()
{
    double total_time = 1000.0;
    double dt = 0.01;

    std::cout << "Running 3D unperturbed orbit tests.\n";
    std::cout << "v_c = " << v_c << ", q = " << q << "\n";

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

        double v_phi = 0.98 * circular_speed(pos_cyl.R);
        double v_R   = 0.02 * v_c;
        double v_z   = 0.02 * v_c;

        Cyl vel_cyl = {v_R, v_phi, v_z};

        Vec3 pos = cyl_to_cart(pos_cyl);
        Vec3 vel = cyl_to_cart_vel(pos_cyl, vel_cyl);

        std::ostringstream basename;
        basename << "orbit3D_" << i;

        simulate_trajectory(pos, vel, basename.str(), total_time, dt);
    }

    std::cout << "Finished all 3D orbit tests.\n";

    return 0;
}