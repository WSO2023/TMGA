

#ifndef CSTCHANGE_GENERATEACHROM_H
#define CSTCHANGE_GENERATEACHROM_H

#include "common.h"
void W_Cal_Average(vector<double>& w);
void C_Cal_Average(vector<vector<double>>& c);
void Calculate_Rank_b(vector<double>& rankList, vector<vector<double>>& c, vector<double>& w);//-w
void IntChr(chromosome& chrom);
vector<int> GnrSS_TS();
vector<int> GnrSS_Lvl();
double GnrMS_Evl(chromosome& chrom);
chromosome GnrChr_HEFT(vector<double> rnk);
chromosome GnrChr_Lvl_EFT();
chromosome GnrChr_TS_EFT();
chromosome GnrChr_TS_Rnd();
void IntPop(vector<vector<chromosome>>& populations, int& stg);
void ChrExc(vector<vector<chromosome>>& populations);
#endif //CSTCHANGE_GENERATEACHROM_H
