#pragma once
#include <vector>
#include <utility>
//#include <iostream>
//#include <fstream>
//using namespace std;

class Arr2d
{
public:
	Arr2d();
	~Arr2d();
	Arr2d(const unsigned int nx, const unsigned int ny);
	/*Arr2d(Arr2d&& arr);*/

   	void Initialize(const unsigned int nx, const unsigned int ny);

	int num_rows();

	int num_cols();

	//void add_data(ofstream& of);

   	double& operator()(const unsigned int x, const unsigned int y);

private:
   	double* arr_ = nullptr;
	unsigned int nx_ = 0;
   	unsigned int ny_ = 0;
};
/*Arr2d arr;
std::vector<Arr2d> allArrays;
*/




