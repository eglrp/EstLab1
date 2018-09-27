
#include "CmnDef.h"
#include "GNSSEKF.h"
#include "MatC.h"
void EKF_Predict(GNSSEKF *& ekf);
void EKF_Execute(GNSSEKF *& ekf);
void EKF_FetchObservation(GNSSEKF *& ekf, double ** satellites_position, int sat_num,
	double * distance,
	double * distance_err);
void EKF_Reset(GNSSEKF *& ekf);
void EKF_FetchResiduals(GNSSEKF *& ekf, double ** satellites_position, int sat_num)
{
	double S[MAX_SATELLITE_NUMBER];
	// 算个残差如何？
	for (int i = 0; i < sat_num; i++)
	{
		S[i] = sqrt(
			pow(ekf->X->data[0][0] - satellites_position[i][0], 2) +
			pow(ekf->X->data[3][0] - satellites_position[i][1], 2) +
			pow(ekf->X->data[6][0] - satellites_position[i][2], 2)
		) + ekf->X->data[9][0];
		ekf->Zp->data[i][0] = S[i];
	}

	mat_minus(ekf->Zp, ekf->Z, ekf->V);
}
//FILE * fp6 = NULL;

GNSSEKF EKFCreate(double DeltaT)
{
	double ttd2 = 0.5 * DeltaT * DeltaT;
	double tttd6 = DeltaT * DeltaT * DeltaT / 6;
	GNSSEKF tot;
	tot.X = malloc_mat(11, 1);
	tot.Xp = malloc_mat(11, 1);
	tot.Dx = eyes(11);
	tot.Dp = eyes(11);

	tot.Z = NULL;
	tot.Dz = NULL;
	tot.H = NULL;
	tot.Ft = NULL;
	tot.Tt = NULL;
	tot.K = NULL;
	tot.V = NULL;
	tot.Zp = NULL;

	//Fai矩阵
	tot.F = malloc_mat(11, 11);
	for (int i = 0; i < 3; i++)
	{
		int offset = i * 3;
		tot.F->data[offset][offset] = 1;
		tot.F->data[offset][offset + 1] = DeltaT;
		tot.F->data[offset][offset + 2] = ttd2;
		tot.F->data[offset + 1][offset + 1] = 1;
		tot.F->data[offset + 1][offset + 2] = DeltaT;
		tot.F->data[offset + 2][offset + 2] = 1;
	}
	tot.F->data[9][9] = 1;
	tot.F->data[9][10] = DeltaT;
	tot.F->data[10][10] = 1;

	//Tao矩阵
	tot.T = malloc_mat(11, 4);
	for (int i = 0; i < 3; i++)
	{
		tot.T->data[3 * i][i] = tttd6;
		tot.T->data[3 * i + 1][i] = ttd2;
		tot.T->data[3 * i + 2][i] = DeltaT;
	}
	tot.T->data[9][3] = ttd2;
	tot.T->data[10][3] = DeltaT;

	tot.De = eyes(4, KF_SYS_NOI);

	mat_trans(tot.F, tot.Ft);
	mat_trans(tot.T, tot.Tt);

	//fp6 = fopen("f6.txt", "w");
	return tot;
}

void EKFProcess(
	GNSSEKF * ekf,
	double * distance,
	double * distance_err,
	double ** satellites_position,
	int sat_num)
{
	free_mat(ekf->V);
	EKF_Predict(ekf);
	EKF_FetchObservation(ekf, satellites_position, sat_num, distance, distance_err);
	EKF_Execute(ekf);
	EKF_FetchResiduals(ekf, satellites_position, sat_num);
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
	EKF_Reset(ekf);
}

void EKF_Execute(GNSSEKF *& ekf)
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
	Matrix * I = eyes(11);
	mat_minus(I, temp2, temp3);
	mat_multiply(temp3, ekf->Dp, ekf->Dx);

	free_mat(temp1); free_mat(temp2); free_mat(temp3);
}
void EKF_Predict(GNSSEKF *& ekf)
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
void EKF_Reset(GNSSEKF *& ekf)
{
	
	free_mat(ekf->K);
}
void EKF_FetchObservation(GNSSEKF *& ekf, double ** satellites_position, int sat_num,
	double * distance,
	double * distance_err)
{
	free_mat(ekf->H);
	free_mat(ekf->Z);
	free_mat(ekf->Dz);
	free_mat(ekf->Zp);

	ekf->H = malloc_mat(sat_num, 11);
	ekf->Z = malloc_mat(sat_num, 1);
	ekf->Dz = malloc_mat(sat_num, sat_num);
	ekf->Zp = malloc_mat(sat_num, 1);

	double x[4] = {
		ekf->Xp->data[0][0],
		ekf->Xp->data[3][0],
		ekf->Xp->data[6][0],
		ekf->Xp->data[9][0]
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
		ekf->Z->data[i][0] = distance[i];
		ekf->Dz->data[i][i] = distance_err[i];

		ekf->H->data[i][0] = -DX0[i] / S[i];
		ekf->H->data[i][3] = -DY0[i] / S[i];
		ekf->H->data[i][6] = -DZ0[i] / S[i];
		ekf->H->data[i][9] = 1;

		ekf->Zp->data[i][0] = S[i] + x[3];
	}
}
