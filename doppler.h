#ifndef DOPPLER_H
#define DOPPLER_H

struct CARandCtx;


/** Compute the number of scattered photons for resonance fluorescence

The number of scattered photons for saturated resonance fluorescence is
given by

nbar = dt S (gamma / (2pi)) (gamma/2)^2 / ((gamma/2)^2 (1 + 2 S) + delta^2)

where S is the laser field intensity in units of the saturation intensity,
gamma is the atomic linewidth, dt is the duration of the time interval,
and delta is the detuning from the atomic transition.
*/

void drift_kick();

double compute_nbar(int i, double inten, double det);

double gaussian_intensity(int i);

double detuning(int i);

double scattering_rate(int i, double inten, double det);



/*
Compute recoil momenta.  This function assumes an isotropic distribution of the
scattered photons.

n        is the number of atoms.
ctx      Random number generator context.
hbar_k   is the momentum of one photon.
nbar     is an array with the number of scattered photons for each atom.
recoils  contains the recoild momenta for each atom in format [px0 py0 pz0 px1
         py1 pz1 ...].
*/
void add_radiation_pressure(int i, struct CARandCtx* ctx, double hbar_k_nrm, double nbar);

static void add_radiation_pressure_one(int i, struct CARandCtx* ctx, double hbar_k_nrm, double nbar);

static void add_radiation_pressure_small_n(int i, struct CARandCtx* ctx, double hbar_k_nrm, int n);

static void add_radiation_pressure_large_n(int i, struct CARandCtx* ctx, double hbar_k_nrm, int n);

void Radiation_Pressure(int i);

#endif

