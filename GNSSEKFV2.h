#pragma once
#pragma once
#include "MatC.h"

struct GNSSEKFV2
{
	Matrix * X;//X Vx Y Vy Z Vz T dT
	Matrix * Xp;
	Matrix * Dx;
	Matrix * Dp;
	Matrix * F;
	Matrix * T;

	Matrix * Ft;
	Matrix * Tt;
	Matrix * De;

	Matrix * Z;
	Matrix * H;
	Matrix * Dz;
	Matrix * Zp;

	Matrix * K;
	Matrix * V;


};
GNSSEKFV2 EKFV2Create(double Dt);
void EKFV2Process(
	GNSSEKFV2 * ekf,
	double * distance,
	double * distance_err,
	double ** satellites_position,
	int sat_num
);

