#ifndef BEND_KICK_UPDATER_H
#define BEND_KICK_UPDATER_H


//void bend_kick(double dt, double Bz, Ensemble ensemble, Force forces, int num_steps);

extern double _epsilon0;
extern double _k;

void ca_bend_kick_update_vector(double time);

void run_trap_potential(double kz);

void run_trap_coulomb(double kz);

void trap_force(double kz);

void coulomb_force();

void rotating_wall_force();

static void coulomb_force_one_pair(int i, const double* r0, const double* r1, double kij);

static double distance(const double *r/*,double delta*/);

void run_no_forces();

void write_to_outfile();

void write_beam_params();

void reset_forces();

#endif
