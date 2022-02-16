#include <cmath>
//doppler cooling of aingle particle in 1D

//Gaussian laser beam parameters
double S0 = 1.0;
double x0 = [1.0,1.0,0.0];
double k0 = [1.0,1.0,0.0];
double sigma = 1.0;
////////

//return intensity at location of particle
double gaussian_intensity(double* x, double S0, double* x0, double* k, double sigma){
	double xp[3];
	double xperp[3];
	double khat[3];
	double kmag = sqrt(k[0]*k[0]+k[1]*k[1]+k[2]*k[2]);
	
	for(int m = 0; m<3; m++){
                khat[i] = k[i]/kmag;
        }
	for(int m = 0; m<3; m++){
		xp[i] = x[i] - x0[i];
	}
	for(int m = 0; m<3; m++){
                xperp[i] = xp[i] - xp[i]*khat[i];
        }
	double xperp_mag = sqrt(xperp[0]*xperp[0]+xperp[1]*xperp[1]+xperp[2]*xperp[2]);

	double inten = S0 * exp(-xperp_mag*xperp_mag/(sigma*sigma));
	return inten;
}

//return detuning of particle from resonance
double detuning(double delta0, double *k, double* v)
{
	double k_dot_v = k[0]*v[0]+k[1]*v[1]+k[2]*v[2];
	return delta0-k_dot_v;
}

void compute_nbar(double dt, double gamma, double inten, double det)
{
	double nbar = dt * scattering_rate(gamma, inten, det);
}

double scattering_rate(double gamma, double inten, double det)
{
	double nu = gamma / (2*M_PI);
	double half_gamma_squared = 0.25*gamma*gamma;
	return inten*nu*half_gamma_squared / (half_gamma_squared*(1.0+2.0*inten)+delta*delta);
}

void add_radiation_pressure(struct CARandCtx* ctx, const double* hbar_k, double nbar)
{ 
	double hbar_k_nrm = sqrt(hbar_k[0]*hbar_k[0]+hbar_k[1]*hbar_k[1]+hbar_k[2]*hbar_k[2]);
	add_radiation_pressure_one(ctx, hbar_k, hbar_k_nrm, nbar);
}

static void add_radiation_pressure_one(struct CARandCtx* ctx, const double* hbar_k, double hbar_k_nrm, double nbar)
{

}
double m = 87*1.67*pow(10,-27); //particle mass
double wavelength = 780.0*pow(10,-9); //laser wavelength
double gamma = 2.0 * M_PI / wavelength;
double hbar = 1.0*pow(10,-34);
double sigma = 1.0*pow(10,-3);
double v = 5;
double dt = 1*

void Radiation_Pressure(double dt, double gamma, double hbar_k, double delta0)
{
	inten = gaussian_intensity(&x, S0, &x0, &k, sigma);
	det = detuning(delta0, &k, &v);
	compute_nbar(dt,gamma,inten,det);
	add_radiation_pressure(rng.context(), hbar_k, nbar);
}
