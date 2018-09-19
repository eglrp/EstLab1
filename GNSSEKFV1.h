#pragma once
#pragma once
#include "MatC.h"

struct GNSSEKFV1
{
	Matrix * X;//X Y Z T
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
GNSSEKFV1 EKFV1Create(double Dt);
void EKFV1Process(
	GNSSEKFV1 * ekf,
	double * distance,
	double * distance_err,
	double ** satellites_position,
	int sat_num
);
