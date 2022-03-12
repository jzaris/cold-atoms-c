#include <cmath>
#include "assert.h"
#include "ca_rand.h"
#include "dSFMT/dSFMT.h"
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "doppler.h"
#include <random>
#include <chrono>
//doppler cooling of aingle particle in 1D

#define SQR(a) ((a) * (a))
#define CA_LARGE_N 128

using namespace std;
using namespace std::chrono;

ofstream outfile;

double x[3] = {0.0,0.0,0.0};
double v[3] = {5.0,0.0,0.0};
double f[3] = {0.0,0.0,0.0};

double m = 87*1.67*pow(10,-27); //particle mass
double wavelength_const = 780.0*pow(10,-9); //laser wavelength
double k_const = 2.0*M_PI/ wavelength_const;
double gamma_const = 2.0 * M_PI * 6.1*pow(10,6);
double hbar = 1.0*pow(10,-34);
double dt = 1*pow(10,-6);

//Gaussian laser beam parameters
const int num_beams = 6;
double x0[num_beams][3] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
double S0[num_beams] = {0.1,0.1,0.1,0.1,0.1,0.1};
double k[num_beams][3] = {k_const,0.0,0.0,-1*k_const,0.0,0.0, 0.0,k_const,0.0,0.0,-1*k_const,0.0, 0.0,0.0,k_const,0.0,0.0,-1*k_const};
double hbar_k[num_beams][3] = {k_const*hbar,0.0,0.0,-1*k_const*hbar,0.0,0.0, 0.0,k_const*hbar,0.0,0.0,-1*k_const*hbar,0.0, 0.0,0.0,k_const*hbar,0.0,0.0,-1*k_const*hbar}; 
double gam[num_beams] = {gamma_const, gamma_const, gamma_const, gamma_const, gamma_const, gamma_const};
double sigma[num_beams] = {1.0*pow(10,-3),1.0*pow(10,-3), 1.0*pow(10,-1),1.0*pow(10,-1), 1.0*pow(10,-1),1.0*pow(10,-1)};  
//double delta0[num_beams] = {-0.5*gam[0],-0.5*gam[1]}; //detuning in Hz
double delta0[num_beams] = {-0.5*gamma_const,-0.5*gamma_const, -0.5*gamma_const,-0.5*gamma_const, -0.5*gamma_const,-0.5*gamma_const}; 
////////

//int seed[num_beams];
int seed;

struct CARandCtx {
        dsfmt_t dsfmt;
};

//struct CARandCtx* ctx[num_beams];
struct CARandCtx* ctx;

int main()
{
	
	std::random_device rd;
        std::mt19937 gen(rd()); //rd is the seed for the mt19937 instance
        std::uniform_int_distribution<std::mt19937::result_type> dist(1,1000000); 
        /*for(int i=0;i<num_beams;i++){
		seed[i] =  dist(gen);
		cout << "seed: " << seed[i];
	}*/
	
	seed = dist(gen);

        /*for (int i =0; i < num_beams; i++){
                ctx[i] = ca_rand_create();
		ca_rand_seed(ctx[i], seed[i]);
		//cout << "seed: " << seed[i];
        }*/
	
	ctx = ca_rand_create();
	ca_rand_seed(ctx, seed);

	
	//cout << k[1][0] << " ";
	//std::random_device rd;
	//std::mt19937 gen(rd());
        //std::normal_distribution<> normal_dist{0.0,1.0};	
	outfile.open("doppler_vel.dat", ios::out | ios::binary);
	double this_vx;
	//outfile << dt << " ";
	auto start = high_resolution_clock::now();
	for(int i=0; i<20000;i++){
		//cout << i;
		//cout << "forcez: " << f[2]; 
		f[0] =0.0;
		f[1]=0.0;
		f[2]=0.0;
		outfile << std::fixed << std::scientific << v[0] << " "; 
	        outfile << std::fixed << std::scientific << v[1] << " ";
	        outfile << std::fixed << std::scientific << v[2] << " ";	
		drift_kick();
		this_vx = v[0];
	}
	cout << "final vx: " << this_vx;
	auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop-start);
        cout << "Time to run simulation: " << duration.count() << endl;
	outfile.close();
	return 1;
}

void drift_kick()
{
	for(int i = 0; i<3; i++)
	{
		x[i] += 0.5*dt*v[i];
	}
	
	for(int i=0; i<num_beams; i++)
	{
		Radiation_Pressure(i);
	}

	for(int i = 0; i<3; i++)
        {
                v[i] += f[i]/m;
        }

	for(int i = 0; i<3; i++)
        {
                x[i] += 0.5*dt*v[i];
        }

}


//return intensity at location of particle
double gaussian_intensity(int i){
	double xp[3];
	double xperp[3];
	double khat[3];
	double kmag = sqrt(k[i][0]*k[i][0]+k[i][1]*k[i][1]+k[i][2]*k[i][2]);
	//cout << kmag << " ";
	for(int m = 0; m<3; m++){
                khat[m] = k[i][m]/kmag;
        }
	//cout << "kx: " << khat[0] << " ky: " << khat[1]  << " kz: " << khat[2] << " ";   
	for(int m = 0; m<3; m++){
		xp[m] = x[m] - x0[i][m];
	}
	for(int m = 0; m<3; m++){
                xperp[m] = xp[m] - abs(xp[m]*khat[m]);
        }
	//cout << xperp[0] << " ";
	double xperp_mag = sqrt(xperp[0]*xperp[0]+xperp[1]*xperp[1]+xperp[2]*xperp[2]);
	//cout << xperp_mag << " ";
	double inten = S0[i] * exp(-xperp_mag*xperp_mag/(sigma[i]*sigma[i]));
	//cout << inten << " ";
	return inten;
}

//return detuning of particle from resonance
double detuning(int i)
{
	double k_dot_v = k[i][0]*v[0]+k[i][1]*v[1]+k[i][2]*v[2];
	//cout << delta0[i] -k_dot_v << " ";
	return delta0[i]-k_dot_v;
}

double compute_nbar(int i, double inten, double det)
{
	//double this_nbar = dt * scattering_rate(i, inten, det);
        //cout << this_nbar << " ";
	//cout << dt << " ";
	//double this_rate = scattering_rate(i, inten,det);
        //cout << this_rate << " ";
	//cout << det<< " ";	
	return dt * scattering_rate(i, inten, det);
}

double scattering_rate(int i, double inten, double det)
{
	double nu = gam[i] / (2*M_PI);
	double half_gamma_squared = 0.25*gam[i]*gam[i];
	//cout << gam[i]<< " ";
	double rate =  inten*nu*half_gamma_squared / (half_gamma_squared*(1.0+2.0*inten)+det*det);
	//cout << "scat. rate: " << rate << " ";
	return rate;
}

void add_radiation_pressure(int i, struct CARandCtx* ctx, double nbar)
{ 
	double hbar_k_nrm = sqrt(hbar_k[i][0]*hbar_k[i][0]+hbar_k[i][1]*hbar_k[i][1]+hbar_k[i][2]*hbar_k[i][2]);
	//cout << hbar_k_nrm << " ";
	add_radiation_pressure_one(i, ctx, hbar_k_nrm, nbar);
}

static void add_radiation_pressure_one(int i, struct CARandCtx* ctx, double hbar_k_nrm, double nbar)
{
	//cout << "add_radiation_pressure_one";
	//cout << nbar << " ";
	int actual_n;
	//std::random_device rd;
	//std::mt19937 gen(rd());
	//std::poisson_distribution<> poisson_dist(nbar);
	

	ca_rand_poisson(ctx, 1, nbar, &actual_n);
	//ca_rand_poisson(ctx, 1, nbar, &actual_n);



	//actual_n = poisson_dist(gen);

	//cout << actual_n;
	//ca_rand_poisson(ctx, 1, nbar, &actual_n);
	//cout << "nbar: " << nbar << " ";
        //ca_rand_poisson(ctx, 1, nbar, &actual_n);
	//cout << actual_n << " ";
        if (actual_n > CA_LARGE_N) {
                add_radiation_pressure_large_n(i, ctx, hbar_k_nrm, actual_n);
        } else {
                add_radiation_pressure_small_n(i, ctx, hbar_k_nrm, actual_n);
        }

}

static void add_radiation_pressure_small_n(int i,
        struct CARandCtx* ctx,
        double hbar_k_nrm,
        int n)
{
        
	//cout << "smalln" << " ";
	//double directions[3][CA_LARGE_N];
	double directions[3][n];
        double nrms[n] = { 0.0 };
        double recoil[3] = { 0.0 };
        int l, j;

	//std::random_device rd;
	//std::mt19937 gen(rd());
	//std::normal_distribution<> normal_dist{0.0,1.0};  
	//cout << n << " ";
        if (0 == n) return;
	//cout << "PASS";
        assert(n <= CA_LARGE_N);
	//cout << "pass ";
	
        ca_rand_gaussian(ctx, n, 0.0, 1.0, &directions[0][0]);
        
	//for(int j = 0; j<n; j++)
	//{
	//	directions[0][j] = normal_dist(gen);
	//	directions[1][j] = normal_dist(gen);
	//	directions[2][j] = normal_dist(gen);
	//}

	ca_rand_gaussian(ctx, n, 0.0, 1.0, &directions[1][0]);
        ca_rand_gaussian(ctx, n, 0.0, 1.0, &directions[2][0]);	
	//ca_rand_gaussian(ctx, n, 0.0, 1.0, &directions[1][0]);
        //ca_rand_gaussian(ctx, n, 0.0, 1.0, &directions[2][0]);
        //ca_rand_gaussian(ctx, n, 0.0, 1.0, &directions[0][0]);  	

        for (l = 0; l < 3; ++l) {
                for (j = 0; j < n; ++j) {
                        nrms[j] += SQR(directions[l][j]);
                }
        }
	//cout << "pass2 ";
        for (j = 0; j < n; ++j) {
                nrms[j] = sqrt(nrms[j]);
        }
	//cout << "pass3 ";
        for (l = 0; l < 3; ++l) {
                for (j = 0; j < n; ++j) {
                        directions[l][j] /= nrms[j];
                }
        }
	//cout << "pass4 ";
        for (l = 0; l < 3; ++l) {
                for (j = 0; j < n; ++j) {
                        recoil[l] += directions[l][j];
                }
                recoil[l] *= hbar_k_nrm;
		//cout << "recx: " << recoil[0] << " recy: " << recoil[1] << " recz: " << recoil[2] << " ";
        }
	//cout << recoil[0] << " ";
	//cout << "recy: " << recoil[1] << " ";
	//cout << "recz: " << recoil[2] << " ";
	//cout << "pass5 ";
        for (l = 0; l < 3; ++l) {
		//cout << n << " " << hbar_k[i][l] << " " << recoil[l] <<" ";
                //cout << n << " ";
		//cout << hbar_k[i][l]<< " ";
		//cout << recoil[l] << " ";
		f[l] += n * hbar_k[i][l] + recoil[l];
        }
}

static void add_radiation_pressure_large_n(int i,
        struct CARandCtx* ctx,
        double hbar_k_nrm,
        int n)
{
        
	//cout << "largen ";
	double recoil[3];
        int l;

        ca_rand_gaussian(ctx, 3, 0.0, hbar_k_nrm * sqrt(n / 3.0), recoil);
        for (l = 0; l < 3; ++l) {
		//cout << n << " " << hbar_k[i][l] << " " << recoil[l] <<" ";
		//cout << n << " ";
                f[l] += n * hbar_k[i][l] + recoil[l];
        }
}



void Radiation_Pressure(int i)
{
	double inten = gaussian_intensity(i); //intensity of ith beam at location of particle
	//cout << inten << " ";
	double det = detuning(i); //detuning of ith beam from atomic resonance, including doppler shift
	double nbar = compute_nbar(i, inten, det);
	//cout << nbar << " ";
	add_radiation_pressure(i, ctx, nbar); //compute momentum kick and add to f (force) array
	//add_radiation_pressure(i, ctx, nbar);
}
