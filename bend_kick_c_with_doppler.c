#include "bend_kick_c_with_doppler.h"
#define GLOBALS_DEF
#include "read_parameters_with_doppler.h"
#undef GLOBALS_DEF
#include "array2d.h"
#include "ca_rand.h"
#include "dSFMT/dSFMT.h"
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <chrono>

#define SQR(a) ((a) * (a))

using namespace std;
using namespace std::chrono;

ofstream outfile;

double _epsilon0 = 8.854187817620 * pow(10,-12);  
double _k = 1.0 / (4.0 * M_PI * _epsilon0);  

int main(int argc, char *argv[]){
	
	//auto start = high_resolution_clock::now();
	initialize(argv[1]); //read input parameters and store values
	//auto stop = high_resolution_clock::now();
	//auto duration = duration_cast<microseconds>(stop-start);
	//cout << "Time to read input and store parameters: " << duration.count() << endl;
	
	/*struct CARandCtx* ctx_test[3];
        for (int i = 0; i< 3;i++)
        {
                //struct CARandCtx* ctxx = ca_rand_create();
                cout << "i: " << i << " ";
                ctx_test[i] = ca_rand_create();
                ca_rand_seed(ctx_test[i], seeds_global(i,0));
        }*/
	cout << "SEED accessed from bend_kick_c_with_doppler.c: " << seed_global << "\n";

	int num_particles = particles_global;
	int num_beams = beams_global;
	//printf("%d ", num_particles);
	//printf("\n");
	//printf("%f ", arr(0,0));
	//printf("%f ", arr(8,0));

	//write scalar parameters to output file
	outfile.open("output.dat", ios::out | ios::binary); //open output file
	outfile << particles_global;
	outfile << " ";
	outfile << std::scientific << steps_global;
	outfile << " ";
	outfile << std::scientific << duration_global;
	outfile << " ";
	outfile << std::scientific <<  dt_global;
	outfile << " ";
	outfile << std::scientific <<  Bz_global;
        outfile << " ";
	outfile << std::scientific <<  kz_global;
        outfile << " ";
	outfile << std::scientific <<  delta_global;
        outfile << " ";
	outfile << std::scientific <<  omegaR_global;
        outfile << " ";
	outfile << std::scientific <<  phi0_global;
        outfile << " ";


	//write array parameters (charge and mass) to output file
	double this_charge;
	for(int i = 0; i<particles_global; i++){
		this_charge = arr(6,i);
		outfile << std::scientific << this_charge;
		outfile << " ";
	}

	double this_mass;
        for(int i = 0; i<particles_global; i++){
                this_mass = arr(7,i);
                outfile << std::scientific << this_mass;
                outfile << " ";
        }

	write_beam_params();
	
	//run_no_forces();
	//run_trap_potential(kz_global);
	auto start = high_resolution_clock::now();
	run_trap_potential(kz_global); //run simulation, writing pos. and vel. to output file at each timestep
	//run_doppler_only();
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop-start);
        cout << "Time to run simulation: " << duration.count() << endl;
	outfile.close(); //close output file
	return 0;
	
}

//run simulation with no forces other than axial field
void run_no_forces(){
        double this_double;

        for(int i = 0; i< steps_global; i++){
                ca_bend_kick_update_vector(dt_global);
                cout << "Particle 1 x-coord: " << arr(0,0) << "\n";
                write_to_outfile();
        }
   
}

//writes current positions and velocities (stored in arr) to output file
void write_to_outfile()
{
	double this_double;
        for(int row = 0; row<6; row++){
        	for(int col = 0; col < particles_global; col++){
                	this_double = arr(row,col);
                        outfile << std::fixed << std::scientific << this_double;
                        outfile << " ";
                }
        }
}

void write_beam_params()
{
        double this_double;
	outfile << beams_global;
        outfile << " ";
        for(int row = 0; row<9; row++){
                for(int col = 0; col < beams_global; col++){
                        this_double = arr_beams(row,col);
                        outfile << std::fixed << std::scientific << this_double;
                        outfile << " ";
                }
        }
}


//sets all force components to zero
void reset_forces()
{
	for(int i = 0; i<3; i++){ 
		for(int j = 0; j<particles_global; j++){
                        forces(i,j) = 0.0;
                }
        }
}

//run simulation with trap force, rot. wall force, and coulomb force
void run_trap_potential(double kz){
	double this_mass;
	
	ca_bend_kick_update_vector(0.5*dt_global); //evolve half timestep
	//cout << "Particle 1 x-coord: " << arr(0,0) << "\n";
        write_to_outfile();

        for(int i = 0; i< steps_global-1; i++){
		trap_force(kz); //update forces array due to momentum kick from trap force in time dt_global
		//cout << "before couloumb: " << forces(2,0) << "\n";
		rotating_wall_force(); //update forces array due to rotating wall force
		coulomb_force(); //update forces array due to coulomb force
		doppler_force();
                //cout << "after couloumb: " << forces(2,0) << "\n";  
		for (int i = 0; i < particles_global; ++i) { //update velocities due to all forces
        		this_mass = arr(7,i);
			arr(3,i) += forces(0,i) / this_mass;  //update x velocity
			arr(4,i) += forces(1,i) / this_mass;  //update y velocity
			arr(5,i) += forces(2,i) / this_mass;  //update z velocity
		}
                ca_bend_kick_update_vector(dt_global);  //update positions and velocities due to axial field
		reset_forces(); //set all forces to zero

                //cout << "Particle 1 x-coord: " << arr(0,0) << "\n";
		write_to_outfile();
        }

	trap_force(kz); //update forces
	rotating_wall_force();
        coulomb_force();	
	doppler_force();
        for (int i = 0; i < particles_global; ++i) { //update velocities
        	this_mass = arr(7,i);
                arr(3,i) += forces(0,i) / this_mass;  
                arr(4,i) += forces(1,i) / this_mass;
                arr(5,i) += forces(2,i) / this_mass;  
	}

	ca_bend_kick_update_vector(0.5*dt_global); //update positions and velocities for last half timestep
        reset_forces();

        //cout << "Particle 1 x-coord: " << arr(0,0) << "\n";
	write_to_outfile();

}

//run simulation with trap force, rot. wall force, and coulomb force
void run_doppler_only(){
        double this_mass;
	
	for(int i = 0; i< particles_global;i++)
	{
		for(int j = 0; j<3;j++)
		{
			arr(j,i)+=0.5*dt_global*arr(j+3,i);
		}	
	}
        write_to_outfile();

        for(int i = 0; i< steps_global-1; i++){
                doppler_force();
                for (int i = 0; i < particles_global; ++i) { //update velocities due to all forces
                        this_mass = arr(7,i);
                        arr(3,i) += forces(0,i) / this_mass;  //update x velocity
                        arr(4,i) += forces(1,i) / this_mass;  //update y velocity
                        arr(5,i) += forces(2,i) / this_mass;  //update z velocity
                }
		for(int i = 0; i< particles_global;i++)
        	{
                	for(int j = 0; j<3;j++)
                	{
                        	arr(j,i)+=0.5*dt_global*arr(j+3,i);
                	}
        	}

                reset_forces(); //set all forces to zero
                write_to_outfile();
        }

        doppler_force();
        for (int i = 0; i < particles_global; ++i) { //update velocities
                this_mass = arr(7,i);
                arr(3,i) += forces(0,i) / this_mass;
                arr(4,i) += forces(1,i) / this_mass;
                arr(5,i) += forces(2,i) / this_mass;
        }
	
	for(int i = 0; i< particles_global;i++)
        {
                for(int j = 0; j<3;j++)
                {
                        arr(j,i)+=0.5*dt_global*arr(j+3,i);
                }
        }

        reset_forces();

        //cout << "Particle 1 x-coord: " << arr(0,0) << "\n";
        write_to_outfile();

}


void ca_bend_kick_update_vector(double time)
{
        double theta, cosTheta, sinTheta, vx_tmp, vy_tmp, omegaB, x, y, z, vx, vy, vz;
        //struct Vec3 *pos = (struct Vec3 *)x;
        //struct Vec3 *vel = (struct Vec3 *)v;

        for (int i = 0; i < particles_global; ++i) {
                omegaB = arr(8,i);
		theta = time * omegaB;
                cosTheta = cos(theta);
                sinTheta = sin(theta);
		vx = arr(3,i);
		vy = arr(4,i);
		vz = arr(5,i);

                arr(0,i) += (sinTheta * vx + (cosTheta - 1.0) * vy) / omegaB;
                arr(1,i) += (-(cosTheta - 1.0) * vx + sinTheta * vy) / omegaB;
                arr(2,i) += time * vz;

                vx_tmp = cosTheta * vx - sinTheta * vy;
                vy_tmp = sinTheta * vx + cosTheta * vy;
                arr(3,i) = vx_tmp;
                arr(4,i) = vy_tmp;
        }
}


void trap_force(double kz)
{
        double alpha, x, y, z;
	double this_charge;

        for (int i = 0; i < particles_global; ++i) {
                this_charge = arr(6,i);
		alpha = this_charge * dt_global;
                x = arr(0,i);
                y = arr(1,i);
                z = arr(2,i);
                forces(0,i) += 0.5 * alpha * kz * x; 
                forces(1,i) += 0.5 * alpha * kz * y;
                forces(2,i) += -alpha * kz * z;
        }

}

void rotating_wall_force()
{
	double alpha, x, y;
        double this_charge;
	double cosphi;
	double sinphi;

	phi_global += 0.5 * omegaR_global *dt_global;
	cosphi = cos(phi_global);
	sinphi = sin(phi_global);

        for (int i = 0; i < particles_global; ++i) {
                this_charge = arr(6,i);
                alpha = this_charge * dt_global;
                x = arr(0,i);
                y = arr(1,i);
                forces(0,i) += alpha*(delta_global*kz_global*(cosphi*cosphi-sinphi*sinphi)*x + 2*cosphi*sinphi*delta_global*kz_global*y);
                forces(1,i) += alpha*(delta_global*kz_global*(sinphi*sinphi-cosphi*cosphi)*y + 2*cosphi*sinphi*delta_global*kz_global*x);
        }
	phi_global += 0.5 * omegaR_global *dt_global;

}


void coulomb_force()
{
	double kp = dt_global * _k;
	//cout << "_k " << _k;
        double ki;
        double r0[3];
	double r1[3];
        double kij;
        int i, j;
	double this_charge, pair_charge;
	
	for (i = 0; i < particles_global; ++i) {
		this_charge = arr(6,i);
                ki = kp * this_charge;
		//cout << "charge " << this_charge << "\n";   
		//cout << "kp " << kp << "\n";   
		//cout << "ki " << ki << "\n";
                r0[0] = arr(0,i); 
		r0[1] = arr(1,i);
		r0[2] = arr(2,i);	
                for (j = 0; j < particles_global; ++j) {
                        if (j == i) continue;
                        r1[0] = arr(0,j);
			r1[1] = arr(1,j);
			r1[2] = arr(2,j);
			pair_charge = arr(6,j);
                        kij = ki * pair_charge;
			//cout << "kij " << kij << "\n";   
                        coulomb_force_one_pair(i,r0,r1,kij);
                }
        }

}

static void coulomb_force_one_pair(int i, const double* r0, const double* r1, double kij)
{
        double r[3];
        double dist, dist_cubed;
        int m;
        for (m = 0; m < 3; ++m) {
                r[m] = r0[m] - r1[m]; //r[0] is distance b/w particles in x-dir, etc.
        } 
        dist = distance(r);
        dist_cubed = dist * dist * dist;
        for (m = 0; m < 3; ++m) {
                forces(m, i) += kij * r[m] / dist_cubed;
		//cout << "kij " << kij << "\n" << "x: " << r[m] << "\n" << "dist: " << dist_cubed << "\n";
        }
}


static double distance(const double *r/*, double delta*/)
{
        double dist = 0.0;
        int v;
        for (v = 0; v < 3; ++v) {
                dist += r[v] * r[v];
        }
        //dist += delta;
        return sqrt(dist);
}






void doppler_force()
{
	//cout << arr(3,0) << " ";
        for(int i = 0; i<particles_global; i++)
        {
                for(int j = 0; j< beams_global; j++)
                {
                        Radiation_Pressure(j,i, dt_global);
                }
        }
}

void Radiation_Pressure(int beam, int ion, double dt)
{
        double inten = gaussian_intensity(beam, ion); //intensity of ith beam at location of particle
        //cout << inten << " ";
        double det = detuning(beam, ion); //detuning of ith beam from atomic resonance, including doppler shift
        double nbar = compute_nbar(dt, inten, det);
	//cout << nbar << " ";
        add_radiation_pressure(beam, ion, ctx_global, nbar); //compute momentum kick and add to f (force) array
        //add_radiation_pressure(i, ctx, nbar);
}

//return intensity at location of particle
double gaussian_intensity(int beam, int ion){
        double xp[3];
        double xperp[3];
        double khat[3];
        double kmag = sqrt(arr_beams(1,beam)*arr_beams(1,beam)+arr_beams(2,beam)*arr_beams(2,beam)+arr_beams(3,beam)*arr_beams(3,beam));
        //cout << kmag << " ";
        for(int m = 0; m<3; m++){
                khat[m] = arr_beams(m+1,beam)/kmag;
        }
        //cout << "kx: " << khat[0] << " ky: " << khat[1]  << " kz: " << khat[2] << " ";
        for(int m = 0; m<3; m++){
                xp[m] = arr(m, ion) - arr_beams(m+5, beam); //ion position - beam displacement
        }
        for(int m = 0; m<3; m++){
                xperp[m] = xp[m] - abs(xp[m]*khat[m]);
        }
        //cout << xperp[0] << " ";
        double xperp_mag = sqrt(xperp[0]*xperp[0]+xperp[1]*xperp[1]+xperp[2]*xperp[2]);
        //cout << xperp_mag << " ";
        double inten = arr_beams(0,beam) * exp(-xperp_mag*xperp_mag/(arr_beams(4,beam)*arr_beams(4,beam)));
        //cout << inten << " ";
        return inten;
}

//return detuning of particle from resonance
double detuning(int beam, int ion)
{
        double k_dot_v = arr_beams(1,beam)*arr(3,ion)+arr_beams(2,beam)*arr(4,ion)+arr_beams(3,beam)*arr(5,ion);
        //cout << arr_beams(8,beam) -k_dot_v << " ";
        return arr_beams(8, beam)-k_dot_v;
}

double compute_nbar(double dt, double inten, double det)
{	
	//double this_nbar = dt * scattering_rate(inten, det);
	//cout << this_nbar << " ";
	//cout << dt << " ";
	//double this_rate = scattering_rate(inten,det);
	//cout << this_rate << " ";
	//cout << det<< " ";
        return dt * scattering_rate(inten, det);
}


double scattering_rate(double inten, double det)
{
        double nu = gamma0_global / (2*M_PI);
	double half_gamma_squared = 0.25*gamma0_global*gamma0_global;
        //cout << gamma0_global<< " ";
	double rate =  inten*nu*half_gamma_squared / (half_gamma_squared*(1.0+2.0*inten)+det*det);
        //cout << "scat. rate: " << rate << " ";
        return rate;
}

void add_radiation_pressure(int beam, int ion, struct CARandCtx* ctx, double nbar)
{
        double hbar_k_nrm = sqrt(hbar_global*hbar_global*(arr_beams(1,beam)*arr_beams(1,beam)+arr_beams(2,beam)*arr_beams(2,beam)+arr_beams(3,beam)*arr_beams(3,beam)));
        //cout << hbar_k_nrm << " ";
        add_radiation_pressure_one(beam, ion,  ctx, hbar_k_nrm, nbar);
}

static void add_radiation_pressure_one(int beam, int ion, struct CARandCtx* ctx, double hbar_k_nrm, double nbar)
{
        int actual_n;
        ca_rand_poisson(ctx, 1, nbar, &actual_n);
        add_radiation_pressure_small_n(beam, ion, ctx, hbar_k_nrm, actual_n);
        

}

static void add_radiation_pressure_small_n(int beam, int ion,
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
	//cout<< n<< " ";
        //std::random_device rd;
        //std::mt19937 gen(rd());
        //std::normal_distribution<> normal_dist{0.0,1.0};
        //cout << n;
        if (0 == n) return;
        //cout << "PASS";
        //assert(n <= CA_LARGE_N);
        //cout << "pass ";

        ca_rand_gaussian(ctx, n, 0.0, 1.0, &directions[0][0]);

        //for(int j = 0; j<n; j++)
        //{
        //      directions[0][j] = normal_dist(gen);
        //      directions[1][j] = normal_dist(gen);
        //      directions[2][j] = normal_dist(gen);
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
                //cout << recoil[0]<< " ";
		//cout << "recx: " << recoil[0] << " recy: " << recoil[1] << " recz: " << recoil[2] << " ";
        }
        //cout << recoil[0] << " ";
        //cout << "recy: " << recoil[1] << " ";
        //cout << "recz: " << recoil[2] << " ";
        //cout << "pass5 ";
        for (l = 0; l < 3; ++l) {
                //cout << n << " " << hbar_k[i][l] << " " << recoil[l] <<" ";
		//cout<< hbar_global*arr_beams(l+1,beam)<< " ";
		//cout << recoil[l] << " ";
                forces(l,ion) += n * hbar_global*arr_beams(l+1, beam) + recoil[l];
        }
}






