#include <cmath>
//doppler cooling simulation in 1D

//Gaussian laser beam parameters
double S0 = 1.0;
double x0 = [1.0,1.0,0.0];
double k0 = [1.0,1.0,0.0];
double sigma = 1.0
////////

double intensity(double* x, double S0, double* x0, double* k, double sigma){
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

	double inten = S0 * exp(-xperp_mag*xperp_mag/(sigma*sigma))
}

double m = 87*1.67*pow(10,-27); //particle mass
double wavelength = 780.0*pow(10,-9); //laser wavelength
double gamma = 2.0 * M_PI / wavelength;
double hbar = 1.0*pow(10,-34);
double sigma = 1.0*pow(10,-3);
double v = 5;
double dt = 1*
