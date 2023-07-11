#include <cstdlib>
#include "GenerateAChrom.h"
#include "GenOperator.h"
#include "tools.hpp"

//{calculate the average execution time of tasks}
void W_Cal_Average(vector<double>& w) {
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        double sum = 0;
        int RscSize = Tasks[i].ElgRsc.size();
        for (int j = 0; j < RscSize; ++j)
            sum += 1.0 / Rscs[Tasks[i].ElgRsc[j]].pc;
        w[i] = Tasks[i].length * sum / RscSize;
    }
}

//{calculate the average transfer time among tasks}
void C_Cal_Average(vector<vector<double>>& c) {
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        if(Tasks[i].parents.size() == 0){
            continue;
        }
        for (int j = 0; j < Tasks[i].parents.size(); ++j) {
            int parent = Tasks[i].parents[j];
            double sum1 = 0;
            double sum = 0;
            sum = ParChildTranFileSizeSum[parent][i] / VALUE;
            for (int k = 0; k < Tasks[i].ElgRsc.size(); ++k) {
                for (int y = 0; y < Tasks[parent].ElgRsc.size(); ++y) {
                    if (Tasks[i].ElgRsc[k] == Tasks[parent].ElgRsc[y]) {
                        continue;
                    } else {
                        sum1 += sum * 8 / XY_MIN(Rscs[Tasks[i].ElgRsc[k]].bw, Rscs[Tasks[parent].ElgRsc[y]].bw);
                    }
                }
            }
            c[parent][i] = sum1 / (double) (Tasks[i].ElgRsc.size() * Tasks[parent].ElgRsc.size());
        }
    }
}

//calculate the rank of tasks based on independent IO using transfer time C[i][j]
void Calculate_Rank_b(vector<double>& RankList, vector<vector<double>>& c, vector<double>& w){
    for (int i = 0; i < TskLstInLvl[TskLstInLvl.size()-1].size(); ++i) {
        int TaskId=TskLstInLvl[TskLstInLvl.size()-1][i];
        RankList[TaskId] = w[TaskId];
    }
    for(int i =TskLstInLvl.size()-2 ;i >=0 ;--i){
        for (int j = 0; j < TskLstInLvl[i].size(); ++j) {
            int TaskId=TskLstInLvl[i][j];
            double ChildMaxRankc = 0;
            for (int k = 0; k < Tasks[TaskId].children.size(); ++k) {
                int tem = Tasks[TaskId].children[k];
                double CompareObject = RankList[tem] + c[TaskId][tem];
                if(ChildMaxRankc  < CompareObject ){
                    ChildMaxRankc = CompareObject;
                }
            }
            RankList[TaskId] =w[TaskId] + ChildMaxRankc;
        }
    }
}

//{initialize chromosome to allocate spaces}
void IntChr(chromosome& chrom) {
    chrom.TskSchLst.resize(comConst.NumOfTsk);
    chrom.RscAlcLst.resize(comConst.NumOfTsk);
    chrom.EndTime.resize(comConst.NumOfTsk);
}

//{generate a topological sort randomly}
vector<int> GnrSS_TS() {
    vector<int> SS;
    vector<int> upr(comConst.NumOfTsk); //the variables for recording the numbers of unscheduled parent tasks
    vector<int> RTI;                    //the set for recording ready tasks whose parent tasks have been scheduled or not exist
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i] = Tasks[i].parents.size();
    }

    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        if (upr[i] == 0) RTI.push_back(i);
    }
    while (!RTI.empty()) {
        int RandVec = rand() % RTI.size();
        int v = RTI[RandVec];
        vector<int>::iterator iter = RTI.begin() + RandVec;
        RTI.erase(iter);
        for (int i = 0; i < Tasks[v].children.size(); ++i) {
            --upr[Tasks[v].children[i]];
            if (upr[Tasks[v].children[i]] == 0) RTI.push_back(Tasks[v].children[i]);
        }
        SS.push_back(v);
    }
    return SS;
}

void IntPop(vector<vector<chromosome>>& Pops,int& stg){
    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<double> ww(comConst.NumOfTsk, 0);
    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk,0));
    W_Cal_Average(ww);                                //calculate the average execution time of tasks
    C_Cal_Average(cc);                                //calculate the average transfer time among tasks
    Calculate_Rank_b(Rank_b,cc, ww);          //calculate the rank of tasks -w
    chromosome chrom_HEFT = GnrChr_HEFT(Rank_b);         //generate a chromosome according to HEFT
    for (int m = 0; m < Parameter_TMGA.NumOfSubPop; ++m) {
        set<chromosome> TemSubPop;                       //use data structure "set" for removing the same chromosome and sorting(temSubPop)
        TemSubPop.insert(chrom_HEFT);
        int nmb = 0;                                     //record the number of attempts
        while(TemSubPop.size() < Parameter_TMGA.NumOfChrInSubPop) {
            ++nmb;
            if(nmb < 2 * Parameter_TMGA.NumOfChrInSubPop) {
                chromosome TemChrom = GnrChr_Lvl_EFT();  //generate a chromosome whose SS is generated randomly based on level and MS is generated based on EFT
                TemSubPop.insert(TemChrom);
            } else if (nmb < 4 * Parameter_TMGA.NumOfChrInSubPop){//-w
                chromosome TemChrom = GnrChr_TS_EFT();   //generate a chromosome whose SS is generated randomly based on TS and MS is generated based on EFT
                TemSubPop.insert(TemChrom);
            } else {
                chromosome TemChrom = GnrChr_TS_Rnd();   //generate a chromosome whose SS is generated random-ly based on TS and MS is generated randomly
                TemSubPop.insert(TemChrom);
            }
            if( nmb >= 4 * Parameter_TMGA.NumOfChrInSubPop){
                stg = 2;
            }
        }
        vector<chromosome> NewSubPop;
        NewSubPop.assign(TemSubPop.begin(),TemSubPop.end());
        Pops.push_back(NewSubPop);
    }
}

void ChrExc(vector<vector<chromosome>>& Pops){
    //{form the elite populaiton EltPop }
    set<chromosome> EltPop;
    for (int m = 0; m < Parameter_TMGA.NumOfSubPop; ++m) {
        EltPop.insert(Pops[m].begin(),Pops[m].end());
    }
    set<chromosome>::iterator iter = EltPop.begin();
    advance(iter,Parameter_TMGA.NumOfEliteOfPop);
    EltPop.erase(iter,EltPop.end());
    //{select the best N different chromosomes form the subpopulaiton and EltPop to form new subpopulation  }
    for (int m = 0; m < Parameter_TMGA.NumOfSubPop; ++m) {
        set<chromosome> NewSubPop ;
        NewSubPop.insert(Pops[m].begin(),Pops[m].end());
        NewSubPop.insert(EltPop.begin(),EltPop.end());
        set<chromosome>::iterator it = NewSubPop.begin();
        advance(it,Parameter_TMGA.NumOfChrInSubPop);
        populations[m].assign(NewSubPop.begin(),it);
    }
}


//{generate a task scheduling order by the levels of tasks from small to large} -xy2
//{Those haveing the same level are ranked arbitrarily among them} -xy2
vector<int> GnrSS_Lvl() {
    vector<int> ch;
    vector<vector<int>> tem = TskLstInLvl;
    for (int i = 0; i < TskLstInLvl.size(); ++i) {
        random_shuffle(tem[i].begin(), tem[i].end());   //arrange the tasks in each level
        for (int j = 0; j < tem[i].size(); ++j) {
            ch.push_back(tem[i][j]);
        }
    }
    return ch;
}

double GnrMS_Evl(chromosome& ch) {
    for (int i = 0; i < comConst.NumOfTsk; ++i)
        ch.RscAlcLst[i] = -1;
    vector<set<double> > ITL;                           //the idle time-slot lists  for all resources
    double makespan = 0;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0);
        a.insert(99999999 * 1.0);
        ITL.push_back(a);
    }
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = -1;
        int TaskIndex = ch.TskSchLst[i];
        double FinalEndTime = 100000000000;
        double FinalStartTime = 0;
        for (int j = 0; j < Tasks[TaskIndex].ElgRsc.size(); ++j) {
            double ReadyTime = 0;
            int v = Tasks[TaskIndex].ElgRsc[j];
            if(Tasks[TaskIndex].parents.size() != 0){
                for (int n = 0; n < Tasks[TaskIndex].parents.size(); ++n) {
                    int ParentIndex = Tasks[TaskIndex].parents[n];
                    int ParentRscIndex = ch.RscAlcLst[ParentIndex];
                    double max = ch.EndTime[ParentIndex];
                    if(v != ParentRscIndex){
                        double TransferData = ParChildTranFileSizeSum[ParentIndex][TaskIndex];
                        max += TransferData / VALUE * 8 / (XY_MIN(Rscs[v].bw,Rscs[ParentRscIndex].bw));
                    }
                    if (ReadyTime < max){
                        ReadyTime = max;
                    }
                }
            }
            double ExeTime = Tasks[TaskIndex].length / Rscs[v].pc;
            double StartTime = 0;
            double EndTime = 0;
            //{Find an idle time-slot as early as possible from ITL}
            set<double>::iterator pre  = ITL[v].begin();
            set<double>::iterator post = ITL[v].begin();
            ++post;
            while(post != ITL[v].end()) {
                if((*post - *pre) >= ExeTime && ReadyTime <= (*post)-ExeTime) {
                    StartTime = XY_MAX(*pre, ReadyTime);
                    break;
                } else {
                    ++pre;
                    ++pre;
                    ++post;
                    ++post;
                }
            }
            EndTime = StartTime + ExeTime;
            //{find/record the earliest finish time}
            if (EndTime < FinalEndTime) {
                FinalStartTime = StartTime;
                FinalEndTime = EndTime;
                RscIndex = v;
            }
        }
        ch.EndTime[TaskIndex] = FinalEndTime;
        ch.RscAlcLst[TaskIndex] = RscIndex;
        //{update ITL}
        if(ITL[RscIndex].find(FinalStartTime) != ITL[RscIndex].end()) {
            ITL[RscIndex].erase(FinalStartTime);
        } else {
            ITL[RscIndex].insert(FinalStartTime);
        }
        if(ITL[RscIndex].find(ch.EndTime[TaskIndex]) != ITL[RscIndex].end()) {
            ITL[RscIndex].erase(ch.EndTime[TaskIndex]);
        } else {
            ITL[RscIndex].insert(ch.EndTime[TaskIndex]);
        }
        makespan = XY_MAX(makespan, FinalEndTime);
    }
    ch.FitnessValue = makespan;
    return makespan;
}

chromosome GnrChr_HEFT(vector<double> Rank_b) {
    vector<int> ind(comConst.NumOfTsk);
    chromosome TemChrom;
    IntChr(TemChrom);
    IndexSort(ind, Rank_b);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        TemChrom.TskSchLst[i] = ind[comConst.NumOfTsk - i - 1];
    }
    GnrMS_Evl(TemChrom);
    return TemChrom;
}

chromosome GnrChr_Lvl_EFT(){
    chromosome TemChrom;
    IntChr(TemChrom);
    TemChrom.TskSchLst = GnrSS_Lvl();
    GnrMS_Evl(TemChrom);
    return TemChrom;
}

chromosome GnrChr_TS_EFT(){
    chromosome TemChrom;
    IntChr(TemChrom);
    TemChrom.TskSchLst = GnrSS_TS();
    GnrMS_Evl(TemChrom);
    return TemChrom;
}

chromosome GnrChr_TS_Rnd(){
    chromosome TemChrom;
    IntChr(TemChrom);
    TemChrom.TskSchLst = GnrSS_TS();
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RandInt = rand() % Tasks[i].ElgRsc.size();
        TemChrom.RscAlcLst[i] = Tasks[i].ElgRsc[RandInt];
    }
    DcdEvl(TemChrom,true);
    return TemChrom;
}
