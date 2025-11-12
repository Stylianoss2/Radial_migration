#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
#include "VecUtils.h"
#include "LogPot.h"

using std::vector;
using std::ofstream;
using std::string;

void simulate_trajectory(Vec3 &orbit_pos, Vec3 &orbit_vel, const double &q, const std::string &basename,
                         double &total_time, double &dt, bool use_pert)
{     
    std::ofstream trajectory_file(basename + "_position.dat");
    std::ofstream parameters_file(basename + "_parameters.dat");
    std::ofstream accel_file(basename + "_accelerations.dat");
    std::ofstream sos_file(basename + "_sos.dat");
    std::ofstream rot_file(basename + "_position_rot.dat");
    
    if (!trajectory_file || !parameters_file || !accel_file || !sos_file || !rot_file) 
    {
        std::cerr << "Error: Cannot open one of the output files.\n";
        return;
    }
    trajectory_file << "# t\tx\ty\tz\tR\tphi\tz\n";
    parameters_file << "# t\tL_z\tKE\tPE\tE\tE_J\n";
    accel_file << "# t\tax_log\tay_log\taz_log\tax_pert\tay_pert\taz_pert\n";
    sos_file << "# t\tpos_x\tvel_x\tpos_y\tvel_y\n"; 
    rot_file << "# t\tx_rot\ty_rot\tz_rot\n";

    int steps = static_cast<int>(total_time / dt);
    int print_interval = steps / 10;

    double prev_z = orbit_pos.z;  //variable to check if orbit crosses equatorial plane

    for (int i = 0; i <= steps; i++)
    {
        double sim_time = static_cast<double>(i) * dt;

        Vec3 log_acc = LogPot_acc(orbit_pos, q);
        Vec3 total_acc = total_acceleration(orbit_pos, q, sim_time);

        double Lz = use_pert ? angular_momentum_inertial_frame(orbit_pos, orbit_vel) : angular_momentum(orbit_pos, orbit_vel);
        double KE = use_pert ? kinetic_energy_inertial_frame(orbit_pos, orbit_vel) : kinetic_energy(orbit_vel);
        double PE = use_pert ? total_potential_energy(orbit_pos, q, sim_time) : potential_energy_unperturbed(orbit_pos, q);
        double E = KE + PE;
        double EJ = use_pert ? jacobi_integral(orbit_pos, orbit_vel, q, sim_time) : E;     


        //////// Find positions in rotating frame//////////////////////////////////////////////////////////////////////
        Cyl pos_cyl = cartesian_to_cylindrical(orbit_pos); //keep position in inertial frame
        Cyl pos_rot = pos_cyl; //Make new variable to store positions in rotational frame
        pos_rot.phi -= omega_p * sim_time;
        if (pos_rot.phi < 0)
        {
            pos_rot.phi += 2.0 * M_PI;
        } 
        if (pos_rot.phi >= 2.0 * M_PI) 
        {
            pos_rot.phi -= 2.0 * M_PI;
        }
        Vec3 pos_rot_cart = cylindrical_to_cartesian(pos_rot);
        rot_file << sim_time << "\t" << pos_rot_cart.x << "\t" << pos_rot_cart.y << "\t" << pos_rot_cart.z << "\n";
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        if (prev_z < 0.0 && orbit_pos.z >= 0.0 && orbit_vel.z > 0.0)
        {
            Cyl cyl_vel = cartesian_to_cylindrical_velocity(orbit_pos, orbit_vel);
            sos_file << sim_time << "\t" << pos_cyl.R << "\t" << cyl_vel.R << "\t" << orbit_pos.y << "\t" << orbit_vel.y << "\t" << "\n"; 
        }
        prev_z = orbit_pos.z;
        
        
        //files
        trajectory_file << sim_time << "\t" << orbit_pos.x << "\t" << orbit_pos.y << "\t" << orbit_pos.z << "\t"
                        << pos_cyl.R << "\t" << pos_cyl.phi << "\t" << pos_cyl.z << "\n";

        parameters_file << sim_time << "\t" << Lz << "\t" << KE << "\t" << PE << "\t" << E << "\t" << EJ << "\t" << "\n";

        accel_file << sim_time << "\t" << log_acc.x << "\t" << log_acc.y << "\t" << log_acc.z << "\t"
                   << total_acc.x << "\t" << total_acc.y << "\t" << total_acc.z << "\n";

        if (i % print_interval == 0)
        {
            std::cout << std::fixed << std::setprecision(5);
            std::cout << "[Step " << i << "] t = " << sim_time << ", L_z = " << Lz << ", E = " << E << ", E_J = " << EJ << "\n";
        }

        if (use_pert)
        {
            Leapfrog_integrator_perturbed(orbit_pos, orbit_vel, q, dt, sim_time);
        }else {
            Leapfrog_integrator_unperturbed(orbit_pos, orbit_vel, q, dt);
        }
    }

    trajectory_file.close();
    parameters_file.close();
    accel_file.close();
    sos_file.close(); 
    rot_file.close();
}

int main(int argc, char* argv[])
{
    bool use_pert = false;

    if (argc > 1)
    {
        std::string arg = argv[1];
        if (arg == "-pert")
        {
            use_pert = true;
        }
        else if (arg == "-no_pert")
        {
            use_pert = false;
        }
        else
        {
        std::cerr << "Unknown option: " << arg << "\n";
        std::cerr << "Write: ./build/bin/LogPotExec -perturbation or -no_perturbation\n";
        return 1;
        }
    }

    double q = 0.7;
    // 1 complete orbit should take 10-100 Myr
    double total_time = 2.0e+3; 
    double dt = 10.0e-4; //Myr

    std::vector<Vec3> orbit_pos = //{ kpc, kpc, kpc}
    {
        {9.0, 3.0, 0.1},
        {5.0, 2.0, 0.1},
        {6.0, 0.0, 0.1}
    };

    std::vector<Vec3> orbit_vel = // { kpc/Myr }
    {
        {0.0, 0.225, 0.0},
        {0.18, 0.20, 0.0},
        {0.0, 0.21, 0.0}
    };

    for (size_t i = 0; i < orbit_pos.size(); i++)
    {
        std::string base = "orbit_" + std::to_string(i + 1);
        simulate_trajectory(orbit_pos[i], orbit_vel[i], q, base, total_time, dt, use_pert);
    }
    return 0;
}

