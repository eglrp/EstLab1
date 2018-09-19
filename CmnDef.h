#pragma once
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#define MAX_SATELLITE_NUMBER 32
#define FREQ1       1.57542E9           /* L1  frequency (Hz) */
#define C           299792458.0         /* speed of light (m/s) */
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
#define Deg (180.0/M_PI)
inline double _fastcall Hopfield(double hgt, double elev)
{
	double t0, p0, e0, h0;
	double t, p, e;
	double dend, elev2, denw, hw, hd, rkd, rkw;
	double trop;

	if (fabs(hgt)>30000.0)   return 0.0;

	t0 = 20 + 273.16;
	p0 = 1013.0;
	e0 = 0.5*exp(-37.2465 + 0.213166*t0 - 0.000256908*t0*t0);
	h0 = 0;
	hw = 11000.0;
	t = t0 - 0.0068*(hgt - h0);
	p = p0*pow(1.0 - 0.0068 / t0*(hgt - h0), 5);
	e = e0*pow((1 - 0.0068 / t0*(hgt - h0)), 2.0)*pow((1.0 - (hgt - h0) / hw), 4.0);
	elev2 = elev*elev * Deg * Deg;
	dend = sqrt(elev2 + 6.25) / Deg;
	denw = sqrt(elev2 + 2.25) / Deg;

	hd = 148.72*t0 - 488.3552;
	rkd = 1.552e-5*p / t*(hd - hgt);
	rkw = 7.46512e-2*(e / t / t)*(hw - hgt);
	trop = (rkd / sin(dend)) + (rkw / sin(denw));
	return trop;
}