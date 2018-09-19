#pragma warning(disable : 4991)

// there are so many things can do with the data.
// 1. observe and compare the Cex from KF & LS
// 2. do MV and RR

#include <stdio.h>
#include <string.h>
#include "GNSSEKF.h"
#include "CmnDef.h"
#include "Valid.h"
#include "MatC.h"

// configurable variables
#define FILE_NAME "satpos_meas.txt"
#define SAVE_NAME "resultkf.txt"
#define SAMPLE_RATE 1 // Hz
#define METHOD "KF"

#define LS_MAX_ITER 20
#define LS_CONV_THRES 0.0001 // m

#define OBS_Sig0 1 // m

// data buffers for independent epochs
satellite sats[MAX_SATELLITE_NUMBER];
double current_solution[4]{ 0,0,0,0 };
double current_time = INVALID;
int satellite_amount = 0;

// saving file
FILE * sf = fopen(SAVE_NAME, "w");

GNSSEKF ekf = EKFCreate(1.0 / SAMPLE_RATE);

bool solve()
{
	if (strcmp(METHOD, "LS") == 0)
	{
		// epoch initialization
		Matrix * Z = malloc_mat(satellite_amount, 1);
		Matrix * H = malloc_mat(satellite_amount, 4);
		Matrix * D = malloc_mat(satellite_amount, satellite_amount);
		Matrix * V = malloc_mat(satellite_amount, 1);
		Matrix * Sig = malloc_mat(4, 4);
		Matrix * X = malloc_mat(4, 1);

		double DX0[MAX_SATELLITE_NUMBER];
		double DY0[MAX_SATELLITE_NUMBER];
		double DZ0[MAX_SATELLITE_NUMBER];
		double S  [MAX_SATELLITE_NUMBER];
		
		double last_solution[4];
		

		for (int i = 0; i < LS_MAX_ITER; i++)
		{
			memcpy(last_solution, current_solution, sizeof(double) * 4);
			for (int j = 0; j < satellite_amount; j++)
			{
				DX0[j] = sats[j].location[0] - current_solution[0];
				DY0[j] = sats[j].location[1] - current_solution[1];
				DZ0[j] = sats[j].location[2] - current_solution[2];
				S[j] = sqrt(DX0[j] * DX0[j] + DY0[j] * DY0[j] + DZ0[j] * DZ0[j]);
			}

			for (int j = 0; j < satellite_amount; j++)
			{
				D->data[j][j] = OBS_Sig0 * OBS_Sig0 / sin(sats[j].elevation);
				Z->data[j][0] = sats[j].pseudorange - S[j] - current_solution[3];
				H->data[j][0] = -DX0[j] / S[j];
				H->data[j][1] = -DY0[j] / S[j];
				H->data[j][2] = -DZ0[j] / S[j];
				H->data[j][3] = 1;
			}

			LMS(Z, H, D, X, Sig, V);

			for(int j = 0; j < 4; j++)
				current_solution[j] += X->data[j][0];

			if (distance(last_solution, current_solution, 3) <= LS_CONV_THRES) {
				// job done
				break;
			}
		}
	}
	else if (strcmp(METHOD, "KF") == 0)
	{
		double distance[MAX_SATELLITE_NUMBER];
		double distance_var[MAX_SATELLITE_NUMBER];
		double ** sat_loc = (double**)alloca(satellite_amount * sizeof(double*));
		for (int i = 0; i < satellite_amount; i++)
		{
			distance[i] = sats[i].pseudorange;
			distance_var[i] = OBS_Sig0 * OBS_Sig0 / sin(sats[i].elevation);

			sat_loc[i] = (double*)alloca(3 * sizeof(double));
			memcpy(sat_loc[i], sats[i].location, sizeof(double) * 3);
		}
		EKFProcess(&ekf, distance, distance_var, sat_loc, satellite_amount);
		current_solution[0] = ekf.X->data[0][0];
		current_solution[1] = ekf.X->data[3][0];
		current_solution[2] = ekf.X->data[6][0];
		current_solution[3] = ekf.X->data[9][0];

		
	}
	fprintf(sf,
		"%lf\t%lf\t%lf\t%lf\n", current_solution[0], current_solution[1], current_solution[2], current_solution[3]
	);
	

	return true;
}

int main()
{

	// open the file
	FILE * fp = fopen(FILE_NAME, "r");
	if (!fp) return INVALID;
	fscanf(fp, "%lf", &current_time);
	if(!is_valid_time(current_time)){
		return INVALID;
	}

	while (!feof(fp))
	{
		int prn = INVALID;
		double time = INVALID;

		satellite & sat = sats[satellite_amount];
		fscanf(fp, "%d%lf%lf%lf%lf%lf%lf\n",
			&sat.prn, &sat.location[0], &sat.location[1], &sat.location[2], &sat.pseudorange, &sat.carrier_phase, &sat.elevation);

		if (feof(fp)) break;

		if (is_valid_sat(&sat)) {
			satellite_amount++;
			sat.elevation /= 180.0 / M_PI; // to arc
		}

		fscanf(fp, "%lf", &time);
		if (!is_valid_time(time)) {
			return INVALID;
		}

		if (current_time != time)
		{
			solve();
			satellite_amount = 0;
			current_time = time;
		}
	}

	_fcloseall();
}