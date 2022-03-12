#ifndef READ_PARAMETERS_H
#define READ_PARAMETERS_H

#include "array2d.h"
#include "ca_rand.h"
#include "dSFMT/dSFMT.h"
#include <vector>
#include <math.h>
//#include <fstream>
//#include <iostream>
//using namespace std;

#ifdef GLOBALS_DEF
#define GLOBALS_EXT 
#else
#define GLOBALS_EXT extern
#endif

/*
struct CARandCtx {
        dsfmt_t dsfmt;
};*/



GLOBALS_EXT double dt_global;
GLOBALS_EXT double duration_global;
GLOBALS_EXT double steps_global;
GLOBALS_EXT int particles_global;
GLOBALS_EXT double Bz_global;
GLOBALS_EXT double kz_global;
GLOBALS_EXT double delta_global;
GLOBALS_EXT double omegaR_global;
GLOBALS_EXT double phi0_global;
GLOBALS_EXT double phi_global;
GLOBALS_EXT Arr2d arr;
GLOBALS_EXT int beams_global;
GLOBALS_EXT Arr2d arr_beams;
GLOBALS_EXT double gamma0_global;
GLOBALS_EXT double hbar_global;
//GLOBALS_EXT Arr2d seeds_global;
GLOBALS_EXT int seed_global;
//GLOBALS_EXT std::vector<CARandCtx*> ctx_global;
GLOBALS_EXT CARandCtx* ctx_global;
GLOBALS_EXT Arr2d forces;
GLOBALS_EXT std::vector<Arr2d> vec;

//GLOBALS_EXT ofstream outfile;



void initialize(const char *filename);

void str_to_double(double *arr, const char *str, int num_particles);

/*#include <vector>

class Arr2d
{
public:
   Arr2d();
   ~Arr2d()
   {
      if (arr)
          delete[] arr;
   }
   Arr2d(const unsigned int nx, const unsigned int ny);
   void Initalize(const unsigned int nx, const unsigned int ny)
   {
      this->nx = nx;
      if (!arr)
      {
	  arr = new double[nx*ny];
      }
   }
   double& operator()(const unsigned int x, const unsigned int y)
   {
      return(arr[nx*y + x]);
   }
private:
   double* m_arr = nullptr;
   unsigned int nx = 0;
   unsigned int ny_ = 0;
};
Arr2d arr;
std::vector<Arr2d> allArrays;

#include "array.h"

Arr2d(const unsigned int nx, const unsigned int ny)
{
   Initialize(nx,ny);
}

void Arr2d::Intiailize(const unsigned int nx, const unsigned int ny)
{
   
}


#include "array.h"

arr.Initialize(10,10);
arr(1,1)++;
allArrays.push_back(arr);
*/

#endif
