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
#define FILE_NAME "satpos_meas4.txt"
#define SAVE_NAME "resultls.txt"
#define ANAL_NAME "ana.txt"
#define LOG_NAME  "log.txt"
#define OBSA_NAME "obs_sta.txt"
#define SAMPLE_RATE 1 // Hz
#define METHOD "LS"

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
FILE * LOG = fopen(LOG_NAME, "w");
FILE * af = fopen(ANAL_NAME, "w");
FILE * of = fopen(OBSA_NAME, "w");

GNSSEKF ekf = EKFCreate(1.0 / SAMPLE_RATE);

double lambda = C / FREQ1;

int ana_prn = 15;
void XYZ2BLH(double * XYZ, double * BLH)
{
	const static double a = 6378137.0;
	const static double F = 1.0 / 298.257223563;
	double e2, Z, dZ, ZdZ, r, sinb, N, x2y2;
	int iter;
	iter = 0;
	r = 0.0;
	N = 0.0;
	sinb = 0.0;
	e2 = 2 * F - F * F;
	x2y2 = XYZ[0] * XYZ[0] + XYZ[1] * XYZ[1];
	dZ = e2 * XYZ[2];
	do
	{
		Z = dZ;
		ZdZ = Z + XYZ[2];
		r = x2y2 + ZdZ*ZdZ;
		sinb = ZdZ / sqrt(r);
		N = a / sqrt(1 - e2*sinb*sinb);
		dZ = N * e2 * sinb;
		iter = iter + 1;
	} while ((iter <= 10) && (fabs(dZ - Z) > 1E-8));
	BLH[0] = atan2(XYZ[1], XYZ[0]);
	BLH[1] = atan2(ZdZ, sqrt(x2y2));
	BLH[2] = sqrt(x2y2 + ZdZ*ZdZ) - N;
}
bool anal()
{
	for (int i = 0; i < satellite_amount; i++)
	{
		if (sats[i].prn == ana_prn) {
			double BLH[3];
			XYZ2BLH(current_solution, BLH);
			satellite & ref = sats[i];
			//if(sats[i].)
			double range = ref.pseudorange - current_solution[3] - Hopfield(BLH[i], ref.elevation);
			double ambi1 = ref.carrier_phase - range / lambda;
			double e_distance = distance(current_solution, ref.location, 3);
			double error_p = range - e_distance;
			double error_l = (ref.carrier_phase - round(ambi1)) * lambda - e_distance;
			fprintf(af, "G%2d %7lf  %12.3lf %12.3lf %13d %16.3lf %12.3lf %12.3lf\n",
				ref.prn, current_time, range, ref.carrier_phase, (int)round(ambi1), ambi1, error_p, error_l);
		}
	}
	return true;
}
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
				fprintf(LOG, "iter: %d, ", i + 1);
				fprintf(LOG, "SPP: %14.3lf, %14.3lf, %14.3lf, %7.3lf, V: ", current_solution[0], current_solution[1], current_solution[2], current_solution[3]);
				for (int i = 0; i < satellite_amount; i++)
					fprintf(LOG, "%6.2lf ", V->data[i][0]);
				
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
	fprintf(LOG, "%7.1lf    ", current_time);
	while (!feof(fp))
	{
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

		fprintf(LOG, "G%2d:%d, ", sat.prn, (int)round(sat.elevation * 180.0 / M_PI / 10));
	
		if (current_time != time)
		{
			fprintf(of, "%d\n", satellite_amount);
			for(int i = 0; i < 14 - satellite_amount;i++)
				fprintf(LOG, "       ");
			fprintf(LOG, "n: %3d, ", satellite_amount);
			if (satellite_amount >= 4)
			{
				solve();
				anal();
			}
			fprintf(LOG, "\n");
			
			satellite_amount = 0;
			current_time = time;
			fprintf(LOG, "%7.1lf    ", current_time);


		}
	}
	for (int i = 0; i < 14 - satellite_amount; i++)
		fprintf(LOG, "       ");
	fprintf(LOG, "n: %3d, ", satellite_amount);
	solve();
	fprintf(LOG, "\n");

	_fcloseall();
}