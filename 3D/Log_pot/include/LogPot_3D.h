#pragma once
#include "VecUtils_3D.h"

// Global potential parameters
extern double v_c;
extern double q;

// Potential and acceleration
double potential_energy(const Vec3& pos);
Vec3 acceleration(const Vec3& pos);

// Energies
double kinetic_energy(const Vec3& vel);
double total_energy(const Vec3& pos, const Vec3& vel);

// Angular momentum
Vec3 angular_momentum(const Vec3& pos, const Vec3& vel);
double angular_momentum_z(const Vec3& pos, const Vec3& vel);

// Circular speed in the midplane
double circular_speed(double R);

// One leapfrog integration step
void leapfrog_step(Vec3& pos, Vec3& vel, double dt);

