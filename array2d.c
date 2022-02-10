#include "array2d.h"
#include "read_parameters.h"
//#include <fstream>
//#include <iostream>

#include <utility>
//using namespace std;

Arr2d::Arr2d()
{
}

Arr2d::~Arr2d()
{
   if (arr_)
      delete[] arr_;
}

/*Arr2d::Arr2d(Arr2d&& arr) :
	nx_(arr.nx_),
	ny_(arr.ny_),
	arr_(std::move(arr.arr_))
{
	arr.arr_ = nullptr;
}*/


Arr2d::Arr2d(const unsigned int nx, const unsigned int ny)
{
	Initialize(nx,ny);
}

void Arr2d::Initialize(const unsigned int nx, const unsigned int ny)
{
	nx_ = nx;
	ny_ = ny;
	if(!arr_){
		arr_ = new double[nx_*ny_];
	}
}

int Arr2d::num_rows(){
	return nx_;
}

int Arr2d::num_cols(){
        return ny_;
}

/*void Arr2d::add_data(ofstream& of)
{
	for(int row = 0; row<9; row++){
		for(int col = 0; col < particles_global; col++){
			of.write(reinterpret_cast<char*>(arr_(row,col)), sizeof(double));
		}
	}
}*/



double& Arr2d::operator()(const unsigned int x, const unsigned int y)
{
	return(arr_[ny_*x+y]);
}	

