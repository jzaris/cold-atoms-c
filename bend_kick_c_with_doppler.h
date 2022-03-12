#ifndef BEND_KICK_UPDATER_H
#define BEND_KICK_UPDATER_H


//void bend_kick(double dt, double Bz, Ensemble ensemble, Force forces, int num_steps);

extern double _epsilon0;
extern double _k;

void ca_bend_kick_update_vector(double time);

void run_trap_potential(double kz);

void run_doppler_only();

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

void doppler_force();

void Radiation_Pressure(int beam, int ion, double dt);

double gaussian_intensity(int beam, int ion);

double detuning(int beam, int ion);

double compute_nbar(double dt, double inten, double det);

double scattering_rate(double inten, double det);

void add_radiation_pressure(int beam, int ion, struct CARandCtx* ctx, double nbar);

static void add_radiation_pressure_one(int beam, int ion, struct CARandCtx* ctx, double hbar_k_nrm, double nbar);

static void add_radiation_pressure_small_n(int beam, int ion,
        struct CARandCtx* ctx,
        double hbar_k_nrm,
        int n);


#endif
