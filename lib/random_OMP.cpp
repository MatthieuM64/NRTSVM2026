/*LIBRARY FOR RANDOM NUMBER GENERATOR*/

#ifndef DEF_RANDOM_MANGEAT_CPP
#define DEF_RANDOM_MANGEAT_CPP

#include <gsl/gsl_randist.h>
#include <cstdint>

const int OMP_MAX_THREADS=omp_get_max_threads();
vector<gsl_rng*> GSL_r(OMP_MAX_THREADS);
vector<gsl_rng*> GSL_r2(OMP_MAX_THREADS);

void init_gsl_ran()
{
	for (int i=0; i<OMP_MAX_THREADS; i++)
	{
		GSL_r[i]=gsl_rng_alloc(gsl_rng_mt19937);
		GSL_r2[i]=gsl_rng_alloc(gsl_rng_mt19937);
	}
}

double ran()
{
	double r=gsl_rng_uniform(GSL_r[omp_get_thread_num()]);
	return (r != 0) ? r : ran();
}

double ran2()
{
	double r=gsl_rng_uniform(GSL_r2[omp_get_thread_num()]);
	return (r != 0) ? r : ran();
}

double gaussian()
{
	static bool eval=false;
	static double next;
	if (not eval)
	{
		double phi=2*M_PI*ran();
		double psi=ran();
		double rad=sqrt(-2*log(psi));
		eval=true;
		next=rad*sin(phi);
		return rad*cos(phi);
	}
	else
	{
		eval=false;
		return next;
	}
}

uint32_t HASH_seed=0;

inline uint32_t fast_hash(uint32_t x)
{
	x ^= x >> 16;
	x *= 0x7feb352d;
	x ^= x >> 15;
	x *= 0x846ca68b;
	x ^= x >> 16;
	return x;
}

inline void set_hash_seed(uint32_t seed)
{
	HASH_seed = seed;
}

inline double ran_hash()
{
	return (fast_hash(HASH_seed++) & 0xFFFFFF) / double(0x1000000);
}

#endif


