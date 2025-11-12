#include <iostream>
#include <cmath>
#include "VecUtils.h"
using std::string;
using std::cos;
using std::sin;
using std::tan;
using std::log;

//Global parameters
inline double m = 2.0;
inline double G = 4.5e-12; // kpc^3 / (Msun * Myr^2)
inline double R_o = 2.5; // kpc
inline double R_CR = 10; // kpc
inline double Sigma_o = 1.226e+9; // Msun / kpc^2
inline double e_s = 0.1;
inline double theta_deg = 15;
inline double theta = theta_deg * M_PI / 180.0;
double alpha = m / std::tan(theta);
inline double omega_p = 0.004; // dφ/dt of pattern  km/s/kpc used 0.004 rad/Myr

inline double v_c = 0.225; //kpc/Myr = 220km/s

//Simple logarithmic potential
Vec3 LogPot_acc(const Vec3 &pos, const double &q) 
{
    double R_squared = pos.x * pos.x + pos.y * pos.y;
    double z_squared = pos.z * pos.z / (q * q);
    double potential = R_squared + z_squared;

    double dPhi_dx = -(v_c * v_c * pos.x) / potential;
    double dPhi_dy = -(v_c * v_c * pos.y) / potential;
    double dPhi_dz = -(v_c * v_c * pos.z) / (q * q * potential);
    
    return {dPhi_dx, dPhi_dy, dPhi_dz};
}

void Leapfrog_integrator_unperturbed(Vec3 &pos, Vec3 &vel, const double &q, const double &dt)
{
    Vec3 r_dash;
    r_dash.x = pos.x + 0.5 * dt * vel.x;
    r_dash.y = pos.y + 0.5 * dt * vel.y;
    r_dash.z = pos.z + 0.5 * dt * vel.z;

    Vec3 accel = LogPot_acc(r_dash, q);
    vel.x += dt * accel.x;
    vel.y += dt * accel.y;
    vel.z += dt * accel.z;

    pos.x = r_dash.x + 0.5 * dt * vel.x;
    pos.y = r_dash.y + 0.5 * dt * vel.y;
    pos.z = r_dash.z + 0.5 * dt * vel.z;
}

double angular_momentum(const Vec3 &pos, const Vec3 &vel)
{
    return  pos.x * vel.y - pos.y * vel.x;
}

double kinetic_energy(const Vec3 &vel)
{
    return 0.5 * (vel.x * vel.x + vel.y * vel.y + vel.z * vel.z);
}

double potential_energy_unperturbed(Vec3 &pos, const double &q)
{
    return 0.5 * v_c * v_c * log(pos.x * pos.x + pos.y * pos.y + (pos.z * pos.z / (q * q)));
}

double total_potential_energy(const Vec3 &cart_pos, const double &q, const double &sim_time)
{
    //The total potential energy in the inertial frame is takenn from Daniel & Wyse (2015) eqn. 36
    Cyl cyl_pos = cartesian_to_cylindrical(cart_pos);
    double log_pot_energy = 0.5 * v_c * v_c * log(cart_pos.x * cart_pos.x + cart_pos.y * cart_pos.y + (cart_pos.z * cart_pos.z / (q * q)));
    double kappa = alpha / cyl_pos.R;
    double Sigma = Sigma_o * exp(-cyl_pos.R/R_o);
    double Phi_S = (2.0 * M_PI * G * e_s * Sigma) / kappa;
    double pert_energy = Phi_S * cos( alpha * log(cyl_pos.R/R_CR) + m*omega_p*sim_time - m*cyl_pos.phi);

    return log_pot_energy + pert_energy;
}

Vec3 total_acceleration(const Vec3 &cart_pos, const double &q, const double &sim_time)
{
    /* 
    Vec3 log_pot_acc; // = LogPot_acc(cart_pos, q);
    Vec3 pert_acc = perturbation_acceleration(cart_pos);
    Vec3 cor_acc = coriolis_acceleration(cart_vel);
    Vec3 cent_acc = centrifugal_acceleration(cart_pos);
    Vec3 total_acc = log_pot_acc + pert_acc; // + cor_acc + cent_acc;
    */
    double delta = 0.000001;
    Vec3 cart_posdx;  Vec3 cart_posdy;  Vec3 cart_posdz;
    Vec3 total_acc;
    Vec3 cart_posdx2; Vec3 cart_posdy2; Vec3 cart_posdz2;
    cart_posdx.x = cart_pos.x + delta;      cart_posdx.y = cart_pos.y;                cart_posdx.z = cart_pos.z;
    cart_posdy.x = cart_pos.x;              cart_posdy.y = cart_pos.y + delta;        cart_posdy.z = cart_pos.z;
    cart_posdz.x = cart_pos.x;              cart_posdz.y = cart_pos.y;                cart_posdz.z = cart_pos.z + delta;
    cart_posdx2.x = cart_pos.x - delta;     cart_posdx2.y = cart_pos.y;               cart_posdx2.z = cart_pos.z;
    cart_posdy2.x = cart_pos.x;             cart_posdy2.y = cart_pos.y - delta;       cart_posdy2.z = cart_pos.z;
    cart_posdz2.x = cart_pos.x;             cart_posdz2.y = cart_pos.y;               cart_posdz2.z = cart_pos.z - delta;
    total_acc.x = (total_potential_energy(cart_posdx2,q,sim_time) - total_potential_energy(cart_posdx,q,sim_time))*0.5/delta;
    total_acc.y = (total_potential_energy(cart_posdy2,q,sim_time) - total_potential_energy(cart_posdy,q,sim_time))*0.5/delta;
    total_acc.z = (total_potential_energy(cart_posdz2,q,sim_time) - total_potential_energy(cart_posdz,q,sim_time))*0.5/delta;
    return total_acc;
}

double effective_potential(const Vec3 &cart_pos, const double &q, const double &sim_time)
{
    //All terms are calculated in the inertial frame of reference
    double phi_total = total_potential_energy(cart_pos, q, sim_time);
    Vec3 omega_p_cart = {0.0, 0.0, omega_p};
    Vec3 cross = cross_product(omega_p_cart, cart_pos);
    double centrifugal_force = 0.5 * (cross.x * cross.x + cross.y * cross.y + cross.z * cross.z);
    return phi_total - centrifugal_force;
}

Vec3 effective_potential_gradient(const Vec3 &cart_pos, const double &q, const double &sim_time)
{
    double delta = 0.000001;
    Vec3 cart_posdx;  Vec3 cart_posdy;  Vec3 cart_posdz;
    Vec3 eff_pot_grad;
    Vec3 cart_posdx2; Vec3 cart_posdy2; Vec3 cart_posdz2;
    cart_posdx.x = cart_pos.x + delta;      cart_posdx.y = cart_pos.y;                cart_posdx.z = cart_pos.z;
    cart_posdy.x = cart_pos.x;              cart_posdy.y = cart_pos.y + delta;        cart_posdy.z = cart_pos.z;
    cart_posdz.x = cart_pos.x;              cart_posdz.y = cart_pos.y;                cart_posdz.z = cart_pos.z + delta;
    cart_posdx2.x = cart_pos.x - delta;     cart_posdx2.y = cart_pos.y;               cart_posdx2.z = cart_pos.z;
    cart_posdy2.x = cart_pos.x;             cart_posdy2.y = cart_pos.y - delta;       cart_posdy2.z = cart_pos.z;
    cart_posdz2.x = cart_pos.x;             cart_posdz2.y = cart_pos.y;               cart_posdz2.z = cart_pos.z - delta;
    eff_pot_grad.x = (effective_potential(cart_posdx2,q,sim_time) - effective_potential(cart_posdx,q,sim_time))*0.5/delta;
    eff_pot_grad.y = (effective_potential(cart_posdy2,q,sim_time) - effective_potential(cart_posdy,q,sim_time))*0.5/delta;
    eff_pot_grad.z = (effective_potential(cart_posdz2,q,sim_time) - effective_potential(cart_posdz,q,sim_time))*0.5/delta;
    return eff_pot_grad;
}

void Leapfrog_integrator_perturbed(Vec3 &pos, Vec3 &vel, const double &q, const double &dt, const double &sim_time)
{
    Vec3 r_dash;
    r_dash.x = pos.x + 0.5 * dt * vel.x;
    r_dash.y = pos.y + 0.5 * dt * vel.y;
    r_dash.z = pos.z + 0.5 * dt * vel.z;

    Vec3 accel = total_acceleration(r_dash, vel, q, sim_time);
    vel.x += dt * accel.x;
    vel.y += dt * accel.y;
    vel.z += dt * accel.z;

    pos.x = r_dash.x + 0.5 * dt * vel.x;
    pos.y = r_dash.y + 0.5 * dt * vel.y;
    pos.z = r_dash.z + 0.5 * dt * vel.z;
}

//Functions to calculate the Jacobi integral in the rotating frame
//Determine the total potential energy in the inertial frame

// velocity in inertial frame of a rotating potential = x_dot + cross_product(pattern speed,cart_pos)
double kinetic_energy_inertial_frame(const Vec3 &cart_pos, const Vec3& cart_vel)
{
    Vec3 omega_p_cart = { 0.0, 0.0, omega_p };
    Vec3 cross_term = cross_product(omega_p_cart, cart_pos);
    Vec3 vel_in = cart_vel + cross_term;

    return 0.5 * (vel_in.x * vel_in.x + vel_in.y * vel_in.y + vel_in.z * vel_in.z);
}

double angular_momentum_inertial_frame(const Vec3 &cart_pos, const Vec3 &cart_vel)
{
    Vec3 omega_p_cart = { 0.0, 0.0, omega_p };
    Vec3 cross_term = {0.0, 0.0, 0.0}; // cross_product(omega_p_cart, cart_pos);
    Vec3 vel_in = cart_vel + cross_term;

    return  cart_pos.x * vel_in.y - cart_pos.y * vel_in.x;
}

double jacobi_integral(const Vec3 &cart_pos, const Vec3& cart_vel, const double &q, const double &sim_time)
{
    double KE = kinetic_energy(cart_vel);
    double PE = total_potential_energy(cart_pos, q, sim_time);
    double Lz = angular_momentum_inertial_frame(cart_pos, cart_vel);

    return KE + PE - (omega_p * Lz);
}














//Rotating frame forces

//Adding a perturbation with constant pattern speed
// Calculate accelerations using formula acceleration = - ∇Φ0 - ∇Φ1 - coriolis - centtrifugal force 
Vec3 coriolis_acceleration(const Vec3& cart_vel)
{
    Vec3 omega_p_cart = { 0.0, 0.0, omega_p};
    Vec3 coriolis_acc = cross_product(omega_p_cart, cart_vel);

    return coriolis_acc * -2.0;
}

Vec3 centrifugal_acceleration(const Vec3 &cart_pos)
{
    Vec3 omega_p_cart = { 0.0, 0.0, omega_p};
    Vec3 intermediate = cross_product(omega_p_cart, cart_pos);
    Vec3 centrifugal_acc = cross_product(omega_p_cart, intermediate);
    return centrifugal_acc * -1.0 ;
}
Vec3 perturbation_acceleration(const Vec3 &cart_pos)
{   
    Cyl cyl_pos = cartesian_to_cylindrical(cart_pos);
    //Determine the potential amplitude and its derivative from Daniel & Wyse (2015) eqn. 37
/*
    double kappa = alpha / cyl_pos.R;
    double Sigma = Sigma_o * exp(-cyl_pos.R/R_o);
    double Phi_S = (2.0 * M_PI * G * e_s * Sigma) / kappa;
    double dPhi_S_dR = (2.0 * M_PI * G * e_s * Sigma * (cyl_pos.R-R_o)) / (R_o * m * (1.0/tan(theta)));

    //Determine partial derivatives of potential perturbation Φ1 in the rotating fram give from Dnaiel & Wyse (2015) eqn. 6
    double dPhi_dR = dPhi_S_dR * cos(m * cyl_pos.phi);
    double dPhi_dphi = -m * Phi_S * sin(m * cyl_pos.phi);
    //Determine derivatives for chain rule
    double dR_dx = cart_pos.x / cyl_pos.R;
    double dR_dy = cart_pos.y / cyl_pos.R;
    double dphi_dx = -cart_pos.y/(cyl_pos.R * cyl_pos.R);
    double dphi_dy = cart_pos.x /(cyl_pos.R * cyl_pos.R);
    //Calculate potential in cartesian coordinates
    double dPhi_dx = (dPhi_dR * dR_dx) + (dPhi_dphi * dphi_dx);
    double dPhi_dy = (dPhi_dR * dR_dy) + (dPhi_dphi * dphi_dy);
    */

    return {0.0, 0.0, 0.0};      //{ -dPhi_dx, -dPhi_dy, 0.0 };
}
