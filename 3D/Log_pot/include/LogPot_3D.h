#pragma once
#include <string>
#include "VecUtils_3D.h"

// Global axisymmetric potential parameters
extern double v_c;
extern double q;

// Global spiral perturbation parameters
extern double m;
extern double G;
extern double R_o;
extern double R_CR;
extern double Sigma_o;
extern double e_s;
extern double theta_deg;
extern double theta;
extern double alpha;
extern double omega_p;
extern double pert_strength;

// Spiral vertical envelope scale height
extern double spiral_z_scale;

// Spiral time envelope parameters
extern double t_spiral_start;
extern double t_spiral_ramp_up;
extern double t_spiral_hold;
extern double t_spiral_ramp_down;

double potential_energy(const Vec3& pos);
Vec3 LogPot_acc(const Vec3& pos);

// Spiral perturbation helpers
double perturbation_scaling_factor(double sim_time);
double spiral_vertical_envelope(double z);
double d_spiral_vertical_envelope_dz(double z);

double spiral_potential_3D(const Vec3& pos, double sim_time);
Vec3 spiral_acceleration_3D(const Vec3& pos, double sim_time);

// Total potential, acceleration, and diagnostics
double total_potential_energy_3D(const Vec3& pos, double sim_time, bool use_pert);
Vec3 total_acceleration_3D(const Vec3& pos, double sim_time, bool use_pert);

double kinetic_energy(const Vec3& vel);
double total_energy(const Vec3& pos, const Vec3& vel);
double total_energy_3D(const Vec3& pos, const Vec3& vel, double sim_time, bool use_pert);
double jacobi_integral_3D(const Vec3& pos, const Vec3& vel, double sim_time, bool use_pert);
double random_energy(const Vec3& pos, const Vec3& vel);
double radial_random_energy(const Vec3& pos, const Vec3& vel);
double vertical_energy(const Vec3& pos, const Vec3& vel);
double radial_action(const Vec3& pos, const Vec3& vel);
double vertical_action(const Vec3& pos, const Vec3& vel);

// Angular momentum
Vec3 angular_momentum(const Vec3& pos, const Vec3& vel);
double angular_momentum_z(const Vec3& pos, const Vec3& vel);

void Leapfrog_integrator_unperturbed(Vec3& pos, Vec3& vel, const double& dt);
void Leapfrog_integrator_3D(Vec3& pos, Vec3& vel, const double& dt, double sim_time, bool use_pert);
void leapfrog_step(Vec3& pos, Vec3& vel, double dt);

// Orbit output
void simulate_trajectory(Vec3& orbit_pos, Vec3& orbit_vel, const std::string& basename, double& total_time, double& dt, bool use_pert);
