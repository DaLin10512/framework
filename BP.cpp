#include"rnbp_a1_a2.h"
#include"const.h"
#include"regidster.h"
#include"bfi.h"
extern int SIZE;

bool is_target(uint64_t x, uint64_t NumTargets, uint64_t *Target);
void PickNewBaseElement_bp(vector<list> &seq, uint64_t NumInputs, uint64_t NumTargets, uint64_t *XorCount, uint64_t * BaseSize, uint64_t *TargetsFound, uint64_t *Target, uint64_t *Base, uint64_t *Dist, uint64_t *NDist, string *Program);

int BP_algorithm(vector<list> &seq, uint64_t NumInputs, uint64_t NumTargets, uint64_t *Target, uint64_t *Dist)
{	
    uint64_t BaseSize = 0;
    uint64_t XorCount = 0;
    uint64_t TargetsFound = 0;
    uint64_t Base[MaxBaseSize] = { 0 };
    uint64_t NDist[MaxBaseSize] = { 0 };

	vector<uint64_t> P;
    vector<uint64_t> Q;
    int randomly_flag = rand() % 2;
    if(randomly_flag == 1)
    {
    	P = generate_p(NumTargets);
    	Q = generate_p(NumInputs);
    	p_m_q(Target, NumTargets, NumInputs, Dist, P, Q);
    }

    string Program[MaxBaseSize];
    InitBase(NumInputs, NumTargets, &BaseSize, &TargetsFound, Base, Dist, Program);
    
    int counter = 0;
    while (TargetsFound < NumTargets)
    {
        counter++;
        if (!EasyMove(seq, NumInputs, NumTargets, &XorCount, &BaseSize, &TargetsFound, Target, Base, Dist, Program)) PickNewBaseElement_bp(seq, NumInputs, NumTargets, &XorCount, &BaseSize, &TargetsFound, Target, Base, Dist, NDist, Program);
     }

	if(randomly_flag == 1)
    {
    	recover_implementation(seq, P, Q);
	}
	
    return XorCount;
}

int TotalDistance_bp(uint64_t NumInputs, uint64_t NumTargets, uint64_t BaseSize, uint64_t NewBase, uint64_t *Target, uint64_t *Base, uint64_t *Dist, uint64_t *NDist)
 { //returns the sum of distances to targets
    int D = 0;
    int t;
    for (unsigned int i = 0; i < NumTargets; i++) 
    {
        t = NewDistance(i, NumInputs, BaseSize, NewBase, Target, Base, Dist);
        NDist[i] = t;
        D = D + t;
    }

    return D;
}

// PickNewBaseElement is only called when there are no 1's in Dist[]
void PickNewBaseElement_bp(vector<list> &seq, uint64_t NumInputs, uint64_t NumTargets, uint64_t *XorCount, uint64_t * BaseSize, uint64_t *TargetsFound, uint64_t *Target, uint64_t *Base, uint64_t *Dist, uint64_t *NDist, string *Program) 
{
    int MinDistance = 0;
    uint64_t TheBest = 0;
    int ThisDist = 0;
    int ThisNorm = 0, OldNorm = 0;
    int besti = 0, bestj = 0, d = 0;
    bool easytarget = false;
    int BestDist[MaxBaseSize] = { 0 };
    uint64_t NewBase = 0;

    MinDistance = (*BaseSize) * NumTargets; //i.e. something big
    OldNorm = 0; //i.e. something small
    //try all pairs of bases
    for (unsigned int i = 0; i < *BaseSize - 1; i++) 
    {
        for (unsigned int j = i+1; j < *BaseSize; j++) 
        {
            NewBase = Base[i] ^ Base[j];
            //sanity check
            if (NewBase == 0) { cout << "a base is 0, should't happen " << endl; exit(0); }
            //if NewBase is not new continue
            if (is_base(NewBase, *BaseSize, Base)) continue;
            //if NewBase is target then choose it
            easytarget = false;
            if (is_target(NewBase, NumTargets, Target))
             {  cout<< "newbase" << NewBase << endl;
                cout << "shouldn't find an easy target here " << endl;
                exit(0);
                easytarget = true;
                besti = i;
                bestj = j;
                TheBest = NewBase;
                break;
            }
            ThisDist = TotalDistance_bp(NumInputs, NumTargets, *BaseSize, NewBase, Target, Base, Dist, NDist); //this also calculates NDist[]
            if (ThisDist <= MinDistance) 
            {
                //calculate Norm
                ThisNorm = 0;
                for (unsigned int k = 0; k < NumTargets; k++) 
                {
                    d = NDist[k];
                    ThisNorm = ThisNorm + d * d;
                }
                //resolve tie in favor of largest norm
                if ((ThisDist < MinDistance) || (ThisNorm > OldNorm)) 
                {
                    besti = i;
                    bestj = j;
                    TheBest = NewBase;
                    for (unsigned int uu = 0; uu < NumTargets; uu++) 
                    {
                        BestDist[uu] = NDist[uu];
                    }
                    MinDistance = ThisDist;
                    OldNorm = ThisNorm;
                }
            }
        }
        if (easytarget) 
        {
            break;
        }
    }
    //update Dist array
    NewBase = TheBest;
    for (unsigned int i = 0; i < NumTargets; i++) 
    {
        Dist[i] = BestDist[i];
    }
    //update Base with TheBest
    Base[*BaseSize] = TheBest;
    //update linear program

    int n = *XorCount + NumTargets; 
    if(NumTargets < NumInputs){ n = *XorCount + NumInputs; }

    string a = Program[besti].substr(0, Program[besti].find(" "));
    string b = Program[bestj].substr(0, Program[bestj].find(" "));
   // Program[BaseSize] = "t" + to_string(XorCount+NumInputs) + " = " + a + " + " + b;
    Program[*BaseSize] = to_string(n) + " = " + a + " + " + b;
	stringstream ss_a, ss_b;
    list row; 
    row.value.push_back(n);
    row.flag = -1;
    row.usd = 0;
    ss_a << a;
    ss_a >> n;
    row.value.push_back(n);
    ss_b << b;
    ss_b >> n;
    row.value.push_back(n);
    seq.push_back(row);

    *BaseSize = *BaseSize + 1;
    *XorCount = *XorCount + 1;
    if (is_target(TheBest, NumTargets, Target)) 
    {
        *TargetsFound = *TargetsFound + 1; //this shouldn't happen
    }
}

bool is_target(uint64_t x, uint64_t NumTargets, uint64_t *Target)
{
    for (unsigned int i = 0; i < NumTargets; i++) 
    {
        if (x == Target[i]) 
        {
            return true;
        }
    }

    return false;
}


