#pragma once

class Rng{

public: 
	Rng();
	~Rng();
	uintptr_t Rng::context();

private:
	CARandCtx* generator_ = nullptr;
}
