#include "Rng.h"
#include "ca_rand.h"
#include "dSFMT/dSFMT.h"

Rng::Rng()
{	
	generator_ = ca_rand_create()
}

Rng::~Rng()
{
   if (generator_)
      delete[] generator_;
}

uintptr_t Rng::context()
{
	return (uintptr_t)generator_;
}


