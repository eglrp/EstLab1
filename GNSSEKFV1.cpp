﻿


#include "CmnDef.h"
#include "GNSSEKFV1.h"
#include "MatC.h"



void EKFV1_Predict(GNSSEKFV1 *& ekf);
void EKFV1_Execute(GNSSEKFV1 *& ekf);
void EKFV1_FetchResiduals(GNSSEKFV1 *& ekf, double ** satellites_position, int sat_num)
{
	double S[MAX_SATELLITE_NUMBER];
	// 算个残差如何？
	for (int i = 0; i < sat_num; i++)
	{
		S[i] = sqrt(
			pow(ekf->X->data[0][0] - satellites_position[i][0], 2) +
			pow(ekf->X->data[1][0] - satellites_position[i][1], 2) +
			pow(ekf->X->data[2][0] - satellites_position[i][2], 2)
		) + ekf->X->data[3][0];
		ekf->Zp->data[i][0] = S[i];
	}

	mat_minus(ekf->Zp, ekf->Z, ekf->V);
}
void EKFV1_FetchObservation(GNSSEKFV1 *& ekf, double ** satellites_position, int sat_num,
	double * distance,
	double * distance_err);
void EKFV1_Reset(GNSSEKFV1 *& ekf);
//FILE * fp6 = NULL;

GNSSEKFV1 EKFV1Create(double DeltaT)
{
	double ttd2 = 0.5 * DeltaT * DeltaT;
	double tttd6 = DeltaT * DeltaT * DeltaT / 6;
	GNSSEKFV1 tot;
	tot.X = malloc_mat(4, 1);
	tot.Xp = malloc_mat(4, 1);
	tot.Dx = eyes(4);
	tot.Dp = eyes(4);

	tot.Z = NULL;
	tot.Dz = NULL;
	tot.H = NULL;
	tot.Ft = NULL;
	tot.Tt = NULL;
	tot.K = NULL;
	tot.V = NULL;
	tot.Zp = NULL;

	//Fai矩阵
	tot.F = eyes(4);

	//Tao矩阵
	tot.T = malloc_mat(4, 4);
	for (int i = 0; i < 3; i++)
	{
		tot.T->data[i][i] = DeltaT;
	}

	tot.De = eyes(4, KF_SYS_NOI);

	mat_trans(tot.F, tot.Ft);
	mat_trans(tot.T, tot.Tt);

	//fp6 = fopen("f6.txt", "w");
	return tot;
}

void EKFV1Process(
	GNSSEKFV1 * ekf,
	double * distance,
	double * distance_err,
	double ** satellites_position,
	int sat_num)
{
	free_mat(ekf->V);
	EKFV1_Predict(ekf);
	EKFV1_FetchObservation(ekf, satellites_position, sat_num, distance, distance_err);
	EKFV1_Execute(ekf);
	EKFV1_FetchResiduals(ekf, satellites_position, sat_num);
	/*fprintf(fp6, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
	ekf->X->data[0][0],
	ekf->X->data[1][0],
	ekf->X->data[2][0],
	ekf->X->data[3][0],
	ekf->X->data[4][0],
	ekf->X->data[5][0],
	ekf->X->data[6][0],
	ekf->X->data[7][0],
	ekf->X->data[8][0],
	ekf->X->data[9][0],
	ekf->X->data[10][0]
	);*/
	EKFV1_Reset(ekf);
}

void EKFV1_Execute(GNSSEKFV1 *& ekf)
{
	//增益矩阵
	Matrix * temp1 = NULL, *temp2 = NULL, *temp3 = NULL, *temp4 = NULL;
	Matrix * Ht = NULL;
	mat_trans(ekf->H, Ht);
	mat_multiply(ekf->Dp, Ht, temp1);
	mat_multiply(ekf->H, ekf->Dp, temp2);
	mat_multiply(temp2, Ht, temp3);
	mat_sum(temp3, ekf->Dz);
	mat_inv(temp3, temp4);
	mat_multiply(temp1, temp4, ekf->K);
	free_mat(temp1); free_mat(temp2); free_mat(temp3); free_mat(temp4);

	//新息序列
	mat_minus(ekf->Z, ekf->Zp, ekf->V);

	//状态滤波
	mat_multiply(ekf->K, ekf->V, temp1);
	mat_sum(ekf->Xp, temp1, ekf->X);

	//滤波方差
	mat_multiply(ekf->K, ekf->H, temp2);
	Matrix * I = eyes(4);
	mat_minus(I, temp2, temp3);
	mat_multiply(temp3, ekf->Dp, ekf->Dx);



	free_mat(temp1); free_mat(temp2); free_mat(temp3);
}
void EKFV1_Predict(GNSSEKFV1 *& ekf)
{
	mat_multiply(ekf->F, ekf->X, ekf->Xp);
	Matrix * temp1 = NULL, *temp2 = NULL, *temp3 = NULL, *temp4 = NULL;
	mat_multiply(ekf->F, ekf->Dx, temp1);
	mat_multiply(temp1, ekf->Ft, temp2);
	mat_multiply(ekf->T, ekf->De, temp3);
	mat_multiply(temp3, ekf->Tt, temp4);
	mat_sum(temp2, temp4, ekf->Dp);

	free_mat(temp1);
	free_mat(temp2);
	free_mat(temp3);
	free_mat(temp4);
}
void EKFV1_Reset(GNSSEKFV1 *& ekf)
{
	//free_mat(ekf->V);
	free_mat(ekf->K);
}
void EKFV1_FetchObservation(GNSSEKFV1 *& ekf, double ** satellites_position, int sat_num,
	double * distance,
	double * distance_err)
{
	free_mat(ekf->H);
	free_mat(ekf->Z);
	free_mat(ekf->Dz);
	free_mat(ekf->Zp);

	ekf->H = malloc_mat(sat_num, 4);
	ekf->Z = malloc_mat(sat_num, 1);
	ekf->Dz = malloc_mat(sat_num, sat_num);
	ekf->Zp = malloc_mat(sat_num, 1);

	double x[4] = {
		ekf->Xp->data[0][0],
		ekf->Xp->data[1][0],
		ekf->Xp->data[2][0],
		ekf->Xp->data[3][0]
	};

	double S[MAX_SATELLITE_NUMBER];
	double DX0[MAX_SATELLITE_NUMBER];
	double DY0[MAX_SATELLITE_NUMBER];
	double DZ0[MAX_SATELLITE_NUMBER];

	for (int i = 0; i < sat_num; i++)
	{
		S[i] = sqrt(
			pow(x[0] - satellites_position[i][0], 2) +
			pow(x[1] - satellites_position[i][1], 2) +
			pow(x[2] - satellites_position[i][2], 2)
		);
		DX0[i] = satellites_position[i][0] - x[0];
		DY0[i] = satellites_position[i][1] - x[1];
		DZ0[i] = satellites_position[i][2] - x[2];
	}
	for (int i = 0; i < sat_num; i++)
	{
		ekf->Z->data[i][0]  = distance[i];
		ekf->Dz->data[i][i] = distance_err[i];

		ekf->H->data[i][0] = -DX0[i] / S[i];
		ekf->H->data[i][1] = -DY0[i] / S[i];
		ekf->H->data[i][2] = -DZ0[i] / S[i];
		ekf->H->data[i][3] = 1;

		ekf->Zp->data[i][0] = S[i] + x[3];
	}
}
