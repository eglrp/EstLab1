#pragma once
#include "CmnDef.h"
#define INVALID -1
bool is_valid_time(double gps_sec)
{
	return (gps_sec >= 0 && gps_sec < 7 * 24 * 3600);
}
//bool is_valid_prn(int prn)
//{
//	return (prn > 0 && prn <= MAX_SATELLITE_NUMBER);
//}
bool is_valid_sat(satellite * sat)
{
	if (sat->elevation <= 0 || sat->elevation > 90)return false;
	if (sat->prn <= 0 || sat->prn > MAX_SATELLITE_NUMBER)return false;
	else return true;
}