#pragma once
#include <string.h>
#include <math.h>
#define MAX_SATELLITE_NUMBER 32
struct satellite {
	int prn;
	double location[3];
	double pseudorange;
	double carrier_phase;
	double elevation;
	satellite()
	{
		memset(location, 0, sizeof(double) * 3);
		pseudorange = 0;
		carrier_phase = 0;
		elevation = 0;
	}
};

inline double _fastcall distance(double * p1, double * p2, int dim = 3)
{
	double tot = 0;
	for (int i = 0; i < dim; i++)
	{
		tot += (p1[i] - p2[i]) * (p1[i] - p2[i]);
	}
	return sqrt(tot);
}