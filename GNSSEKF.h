#pragma once
#include "MatC.h"

struct GNSSEKF
{
	Matrix * X;//X Vx Ax Y Vy Ay Z Vz Az T dT
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
GNSSEKF EKFCreate(double Dt);
void EKFProcess(
	GNSSEKF * ekf,
	double * distance,
	double * distance_err,
	double ** satellites_position,
	int sat_num
);

