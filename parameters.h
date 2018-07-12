#ifndef PARAMETERS_H_INCLUDED
#define PARAMETERS_H_INCLUDED

/**
 * @brief Structures of parameters
 *
 **/

struct Parameters
{
	unsigned tau_2D;
	unsigned tau_3D;
	unsigned N;
	unsigned Nf;
	unsigned Ns;
	unsigned Npr;
	unsigned Nb;
	unsigned k;
	unsigned p;
	float d;
	float lambda3D;
	float tau;

	int n = 16;
};

#endif
