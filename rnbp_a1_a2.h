#ifndef _BP_H_
#define _BP_H_
#include"const.h"
#include"regidster.h"
#include"bfi.h"
//typedef unsigned long long int uint64_t;//for linux
const int MaxBaseSize = 1000;
const bool PRuint64_tROWS = true;
        
#define SIZE2 64
#define LOOP_SIZE 100
#define LARGE 100000
#define computeSum(table,table_size,s) s=0; for(uint64_t sum_i=0; sum_i<table_size; sum_i++) s+=table[sum_i];

struct Element
{
    uint64_t parent_i;
    uint64_t parent_j;
    uint64_t newDist[SIZE2];
};

int A1_A2_RNBP_algorithm(vector<list> &seq, uint64_t NumInputs, uint64_t NumTargets, int no_choice, uint64_t *Target, uint64_t *Dist, uint64_t *InitDist);
int Paar_algorithm1(vector<list> &seq, uint64_t NumInputs, uint64_t NumTargets, uint64_t *Target);
int RPaar_algorithm1(vector<list> &seq, uint64_t NumInputs, uint64_t NumTargets, uint64_t *Target);
int BP_algorithm(vector<list> &seq, uint64_t NumInputs, uint64_t NumTargets, uint64_t *Target, uint64_t *Dist);
void InitBase(uint64_t NumInputs, uint64_t NumTargets, uint64_t *BaseSize, uint64_t *TargetsFound, uint64_t *Base, uint64_t *Dist, string *Program);
bool is_base(uint64_t x, uint64_t BaseSize, uint64_t *Base);
uint64_t NewDistance(uint64_t u, uint64_t NumInputs, uint64_t BaseSize, uint64_t NewBase, uint64_t *Target, uint64_t *Base, uint64_t *Dist); //calculates the distance from the base to Target[u]
void TotalDistance(uint64_t NumTargets, uint64_t NumInputs, uint64_t BaseSize, uint64_t NewBase, uint64_t *Target, uint64_t *Base, uint64_t *Dist, uint64_t *NDist); //returns the sum of distances to targets
bool reachable(uint64_t T, uint64_t K, uint64_t S, uint64_t BaseSize, uint64_t *Base);
bool EasyMove(vector<list> &seq, uint64_t NumInputs, uint64_t NumTargets, uint64_t *XorCount, uint64_t *BaseSize, uint64_t *TargetsFound, uint64_t *Target, uint64_t *Base, uint64_t *Dist, string *Program); //if any two bases add up to a target, pick them
void PickNewBaseElement(vector<list> &seq, int no_alg, uint64_t NumInputs, uint64_t NumTargets, uint64_t *XorCount, uint64_t *BaseSize, uint64_t TargetsFound, uint64_t *Target, uint64_t *Base, uint64_t *Dist, uint64_t *NDist, string *Program);
void refreshDist(uint64_t NumTargets, uint64_t *Dist, uint64_t *InitDist); // refresh the distance to the targets for subsequent rounds
uint64_t RNBP(Element AllElements[], uint64_t counter, uint64_t NumInputs, uint64_t NumTargets); // RNBP
uint64_t A1(Element AllElements[], uint64_t counter, uint64_t NumInputs, uint64_t NumTargets, uint64_t TargetsFound, uint64_t *Dist); // A1
uint64_t A2(Element AllElements[], uint64_t counter, uint64_t NumInputs, uint64_t NumTargets, uint64_t TargetsFound, uint64_t *Dist); // A2
uint64_t calculateDist(uint64_t A[],uint64_t length);
uint64_t calculateNorm(uint64_t A[],uint64_t length);
bool filtering(uint64_t tempDist[], vector<uint64_t> filter_indices, uint64_t *Dist);



#endif




