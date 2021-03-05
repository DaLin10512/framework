#include"rnbp_a1_a2.h"
#include"regidster.h"
extern mt19937 rand_generator;
extern int SIZE;

int A1_A2_RNBP_algorithm(vector<list> &seq, uint64_t NumInputs, uint64_t NumTargets, int no_choice, uint64_t *Target, uint64_t *Dist, uint64_t *InitDist)
{	
	vector<uint64_t> P;
    vector<uint64_t> Q;
    int randomly_flag = rand() % 2;
    if(randomly_flag == 1)
    {
		P = generate_p(NumTargets);
    	Q = generate_p(NumInputs);
    	p_m_q(Target, NumTargets, NumInputs, Dist, P, Q);
    }
  
 	for(int i = 0; i < NumTargets; i++)
    {
        InitDist[i] = Dist[i];
    }

	int BestCount = LARGE;
    //clock_t start = clock();
	//while (((clock()-start)/CLOCKS_PER_SEC < TIME_LIMIT)) 
	int round = 0;
    while(round < 1)
    {	
    	round++;
		vector<list> t_seq;
        uint64_t XorCount = 0;
        uint64_t BaseSize= 0;
        uint64_t TargetsFound = 0;
        uint64_t Base[MaxBaseSize] = { 0 };
        uint64_t NDist[MaxBaseSize] = { 0 };
        string Program[MaxBaseSize];
        // refreshing the distance and base for subsequent rounds
        refreshDist(NumTargets, Dist, InitDist);
        InitBase(NumInputs, NumTargets, &BaseSize, &TargetsFound, Base, Dist, Program);
        // main loop
        while ((TargetsFound < NumTargets))
        {
            if (!EasyMove(t_seq, NumInputs, NumTargets, &XorCount, &BaseSize, &TargetsFound, Target, Base, Dist, Program)) PickNewBaseElement(t_seq, no_choice, NumInputs, NumTargets, &XorCount, &BaseSize, TargetsFound, Target, Base, Dist, NDist, Program);
        }
 
        if (XorCount < BestCount)
        {	
            BestCount = XorCount;
            seq.erase(seq.begin() + NumInputs, seq.end());
            seq.insert(seq.begin() + NumInputs, t_seq.begin(), t_seq.end());
    	}
    }

    if(randomly_flag == 1)
    {
    	recover_implementation(seq, P, Q);
	}
	return BestCount;
}
vector<list> Algorithm_B(vector<list> e_seq, vector<int> input, uint64_t *Target, uint64_t *Dist, uint64_t *InitDist)
{	
	vector<list> tmp_seq;
    int SIZE_B = input.size();
    uint64_t NumInputs = input.size();
    uint64_t NumTargets = e_seq.size();
	
    for(int i = 0; i < SIZE_B; i++)
	{
		list tmp;
		tmp.value.push_back(i);
		tmp.flag = -1;
		tmp.usd = 0;
		tmp_seq.push_back(tmp); 
	 }
	rand_generator.seed(time(0));
	uniform_int_distribution<int> rand_distribution(0, 5);
	int OPTION_B = rand_distribution(rand_generator);

	if(OPTION_B >= 3)
	{
		A1_A2_RNBP_algorithm(tmp_seq, NumInputs, NumTargets, OPTION_B, Target, Dist, InitDist);
	}
	else if(OPTION_B == 0)
	{
		Paar_algorithm1(tmp_seq, NumInputs, NumTargets, Target);
	}
	else if(OPTION_B == 1)
	{
		RPaar_algorithm1(tmp_seq, NumInputs, NumTargets, Target);
	}
    else if(OPTION_B == 2)
    {
        BP_algorithm(tmp_seq, NumInputs, NumTargets, Target, Dist);
    }
    else
    {
        cout << "No choice of algorithm B!" << endl;
        exit(0);
    }

    //cout<<"          equi_seq length here is :  "<<tmp_seq.size()-SIZE_B<<endl;

	int N = 1000;

	for(int i = SIZE_B; i < tmp_seq.size(); i++)
	{//The index of the new variants idx' = 400+ idx; 
		for(int j = i + 1; j < tmp_seq.size(); j++)
		{	
			for(int k = 1; k < tmp_seq[j].value.size(); k++)
			{
				if(tmp_seq[j].value[k] == tmp_seq[i].value[0])
				{
					tmp_seq[j].value[k] = tmp_seq[i].value[0] + N;
				}
			}	
		}
		tmp_seq[i].value[0] = tmp_seq[i].value[0] + N;
    }
    for(int i = SIZE_B; i < tmp_seq.size(); i++)
    {
		for(int j = 1; j < tmp_seq[i].value.size(); j++)
		{
			if(tmp_seq[i].value[j] < NumInputs)
			{
				tmp_seq[i].value[j] = input[tmp_seq[i].value[j]];
			}
		}
	}

	for(int i = SIZE_B; i < tmp_seq.size(); i++)
	{
		if(tmp_seq[i].flag >= 0)
		{
			for(int j = i + 1; j < tmp_seq.size(); j++)
			{
				for(int k = 1; k < tmp_seq[j].value.size(); k++)
				{
					if(tmp_seq[j].value[k] == tmp_seq[i].value[0])
					{
						tmp_seq[j].value[k] = e_seq[tmp_seq[i].flag].value[0];
					}
				}
			}
			tmp_seq[i].value[0] = e_seq[tmp_seq[i].flag].value[0];
			if(e_seq[tmp_seq[i].flag].flag >= 0)
			{
				tmp_seq[i].flag = e_seq[tmp_seq[i].flag].flag;
			}
            else if(e_seq[tmp_seq[i].flag].flag < 0)
            {
                tmp_seq[i].flag = -1;
            }
		}
	}
	return tmp_seq;
}

void InitBase(uint64_t NumInputs, uint64_t NumTargets, uint64_t *BaseSize, uint64_t *TargetsFound, uint64_t *Base, uint64_t *Dist, string *Program) 
{
    *TargetsFound = 0;
    Base[0] = 1;
    //Program[0] = "x0";
    Program[0] = "0";
    stringstream ss;
    string s;
    for (uint64_t i = 1; i < NumInputs; i++) 
    {
        ss << i;
        s = ss.str();
        Base[i] = 2 * Base[i - 1];
        //Program[i] = "x" + s;
        Program[i] = s;
        ss.str("");
    }
    *BaseSize = NumInputs; //initial base is just the xi's
    for (uint64_t i = 0; i < NumTargets; i++)
	{
        if (Dist[i] == 0) 
		{
            *TargetsFound = *TargetsFound + 1;
        }
    }
}

void TotalDistance(uint64_t NumTargets, uint64_t NumInputs, uint64_t BaseSize, uint64_t NewBase, uint64_t *Target, uint64_t *Base, uint64_t *Dist, uint64_t *NDist) 
{ //returns the sum of distances to targets
    uint64_t D = 0;
    uint64_t t;
    for (uint64_t i = 0; i < NumTargets; i++) 
    {
        t = NewDistance(i, NumInputs, BaseSize, NewBase, Target, Base, Dist);
        NDist[i] = t;
    }
}

bool EasyMove(vector<list> &seq, uint64_t NumInputs, uint64_t NumTargets, uint64_t *XorCount, uint64_t *BaseSize, uint64_t *TargetsFound, uint64_t *Target, uint64_t *Base, uint64_t *Dist, string *Program) 
{
    uint64_t t;
    bool foundone = false;
    /*
    for(uint64_t i = 0; i < NumTargets; i++) 
    {
        cout<<"Dist_"<<i<<Dist[i]<<endl;
    }
    */
	//cout << "Entry EasyMove1" << endl;
    //see if anything in the distance vector is 1
    for(uint64_t i = 0; i < NumTargets; i++) 
    {
        if (Dist[i] == 1) 
        {
            foundone = true;
            t = i;
            break;
        }
    }

    if (!foundone)
    {
        return false;
    }
    //update Dist array
	//cout<<"Entry EasyMove2"<<endl;
    uint64_t NewBase = Target[t];
	//cout<<"MunTargets "<<NumTargets<<endl;
    for (int u = 0; u < NumTargets; u++) 
	{
		//cout<<"Entry EasyMove3"<<endl;
        Dist[u] = NewDistance(u, NumInputs, *BaseSize, NewBase, Target, Base, Dist);
    }
    /*
    for(int i = 0; i < *BaseSize; i++)
    {
    	cout << "base_" << i << " = " << Base[i] << endl;
    }
    */
	//cout << "Entry EasyMove3" << endl;
    //update Base with NewBase
    Base[*BaseSize] = NewBase;
    //find which lines in Base caused this
    string a,b;
    for (uint64_t i = 0; i < *BaseSize; i++) 
    {
        for (uint64_t j = i + 1; j < *BaseSize; j++) 
        {
            if ((Base[i] ^ Base[j]) == Target[t]) 
            {	
            	//cout << "find Target[t]" << Target[t] << endl;
                a = Program[i].substr(0, Program[i].find(" "));
                b = Program[j].substr(0, Program[j].find(" "));
                break;
            }
        }
    }
    stringstream ss, ss_a, ss_b, ss_s1;
    string s1;
    int idx = *XorCount + NumTargets;
    if(NumTargets < NumInputs)
    {
    	idx = *XorCount + NumInputs;
    }  
    ss << idx;
    s1 = ss.str();
    ss.str("");
    //Program[BaseSize] = "y" + s1 + " = " + a + " + " + b;
    Program[*BaseSize] = s1 + " = " + a + " + " + b;

    list row;
    int n = NumTargets + seq.size();
    if(NumTargets< NumInputs){n = NumInputs + seq.size();}
    row.value.push_back(idx);
    ss_a << a;
    ss_a >> n;
    row.value.push_back(n);
    ss_b << b;
    ss_b >> n;
    row.value.push_back(n);
    row.flag  = t;
    row.usd = 0;
    seq.push_back(row);

    *BaseSize = *BaseSize + 1;
    *XorCount = *XorCount + 1;
    *TargetsFound = *TargetsFound + 1;
    return true;
}

// PickNewBaseElement is only called when there are no 1's in Dist[]
void PickNewBaseElement(vector<list> &seq, int no_alg, uint64_t NumInputs, uint64_t NumTargets, uint64_t *XorCount, uint64_t *BaseSize, uint64_t TargetsFound, uint64_t *Target, uint64_t *Base, uint64_t *Dist, uint64_t *NDist,  string *Program) 
{
    // Allocate memory for all possible bases
    Element* AllElements = new Element[(*BaseSize) * ((*BaseSize) - 1)];
    uint64_t counter = 0; // counter to track last element

    for (uint64_t i = 0; i < *BaseSize - 1; i++) 
    {
        for (uint64_t j = i+1; j < *BaseSize; j++) 
        {
            uint64_t NewBase = Base[i] ^ Base[j];
            TotalDistance(NumTargets, NumInputs, *BaseSize, NewBase, Target, Base, Dist, NDist); //this calculates NDist[]

            // Putting in the data uint64_to the array
            for (uint64_t k = 0; k < NumTargets; k++)
            {
                AllElements[counter].newDist[k] = NDist[k];
            }
            AllElements[counter].parent_i = i;
            AllElements[counter].parent_j = j;
            counter++;
        }
    }
    uint64_t chosen;
    // selecting the best pair to XOR
    if(no_alg == 3){chosen = RNBP(AllElements, counter, NumInputs, NumTargets);}
    else if(no_alg == 4){chosen = A1(AllElements, counter, NumInputs, NumTargets, TargetsFound, Dist);}
    else if(no_alg == 5){chosen = A2(AllElements, counter, NumInputs, NumTargets, TargetsFound, Dist);}
    else 
    	{
    		cout << "No choice of the algorithm here！" << endl;
			getchar();
		}

    // Update using the result returned by the criteria
    uint64_t bestparent_i = AllElements[chosen].parent_i;
    uint64_t bestparent_j = AllElements[chosen].parent_j;
    Base[*BaseSize] = Base[bestparent_i] ^ Base[bestparent_j];

    // Update the Dist Array
    for (uint64_t i = 0; i < NumTargets; i++)
    {
        Dist[i] = AllElements[chosen].newDist[i];
    }

    string a = Program[bestparent_i].substr(0, Program[bestparent_i].find(" "));
    string b = Program[bestparent_j].substr(0, Program[bestparent_j].find(" "));
    stringstream ss, ss_a, ss_b;
    string s2;
    int idx = *XorCount + NumTargets;
    if(NumTargets < NumInputs){idx = *XorCount + NumInputs;}
    ss << idx;
    s2 = ss.str();
    ss.str("");
    //Program[BaseSize] = "t" + s2 + " = " + a + " + " + b;
    Program[*BaseSize] = s2 + " = " + a + " + " + b;

    int n = *XorCount + NumTargets; 
    if(NumTargets < NumInputs){ n = *XorCount + NumInputs; }
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
    //free up the memory
    free(AllElements);
    return;
}

// Original BP random criteria
uint64_t RNBP(Element AllElements[], uint64_t counter, uint64_t NumInputs, uint64_t NumTargets)
{
    // initialization
    uint64_t bestDist = LARGE;
    uint64_t bestNorm = -1 * LARGE;
    uint64_t currentDist, currentNorm;
    vector<uint64_t> candidates;
    for (uint64_t i = 0; i < counter; i++)
    {
        currentDist = calculateDist(AllElements[i].newDist,NumTargets);
        currentNorm = calculateNorm(AllElements[i].newDist,NumTargets);
        if ((currentDist < bestDist) || (currentDist == bestDist && currentNorm > bestNorm))
        {
            // Updating the best distance and norm
            bestDist = currentDist;
            bestNorm = currentNorm;
            // clear previous candidates
            candidates.clear();
            // inputting the new candidates
            candidates.push_back(i);
        }
        else if (currentDist == bestDist && currentNorm == bestNorm)
        {
            // equal candidates
            candidates.push_back(i);
        }
    }

    // randomly choose one of the candidates
    rand_generator.seed(time(0));
    uniform_int_distribution<int> rand_distribution(0, candidates.size() - 1);
    int rand_num = rand_distribution(rand_generator);
    return candidates[rand_num];
}

uint64_t A1(Element AllElements[], uint64_t counter, uint64_t NumInputs, uint64_t NumTargets, uint64_t TargetsFound, uint64_t *Dist)
{
    // Applying the Filter
    uint64_t nearest = 1; // change this to relax the filter of nearest targets
    uint64_t filter_dist; // keep track of the largest distance that will pass through the filter
    vector<uint64_t> filter_indices; // keep track of the indices of Target that satisfy the filter

    // sort the distance array
    uint64_t sorted_dist[NumTargets - TargetsFound];
    uint64_t next_index = 0;
    for (uint64_t i = 0; i < NumTargets; i++)
    {
        if (Dist[i] == 0) continue;
        sorted_dist[next_index++] = Dist[i];
    }
    sort(sorted_dist, sorted_dist + NumTargets - TargetsFound);

    // indices that can pass through the filter
    filter_dist = sorted_dist[min(nearest - 1, NumTargets - TargetsFound - 1)]; // largest distance that will pass through the filter
    for (uint64_t i = 0; i < NumTargets; i++)
    {

        if (Dist[i] <= filter_dist && Dist[i] > 0) 
        {
            filter_indices.push_back(i);
        }
    }

    // initialization
    uint64_t bestDist = LARGE;
    uint64_t bestNorm = -1 * LARGE;
    uint64_t currentDist, currentNorm;
    vector<uint64_t> candidates;


    for (uint64_t i = 0; i < counter; i++)
    {
        // Filtering
        if (!filtering(AllElements[i].newDist, filter_indices, Dist)) continue;
        // Normal BP rand
        currentDist = calculateDist(AllElements[i].newDist, NumTargets);
        currentNorm = calculateNorm(AllElements[i].newDist, NumTargets);
        if ((currentDist < bestDist) || (currentDist == bestDist && currentNorm > bestNorm))
        {
            // Updating the best distance and norm
            bestDist = currentDist;
            bestNorm = currentNorm;
            // clear previous candidates
            candidates.clear();
            // inputting the new candidates
            candidates.push_back(i);
        }
        else if (currentDist == bestDist && currentNorm == bestNorm)
        {
            // equal candidates
            candidates.push_back(i);
        }
    }

    // randomly choose one of the candidates
    rand_generator.seed(time(0));
    uniform_int_distribution<int> rand_distribution(0, candidates.size() - 1);
    int rand_num = rand_distribution(rand_generator);
    return candidates[rand_num];
}


// BP rand with the filter of choosing the nearest target
uint64_t A2(Element AllElements[], uint64_t counter, uint64_t NumInputs, uint64_t NumTargets, uint64_t TargetsFound, uint64_t *Dist)
{
    // Applying the Filter
    uint64_t nearest = 1; // change this to relax the filter of nearest targets
    uint64_t filter_dist; // keep track of the largest distance that will pass through the filter
    vector<uint64_t> filter_indices; // keep track of the indices of Target that satisfy the filter

    // sort the distance array
    uint64_t sorted_dist[NumTargets - TargetsFound];
    uint64_t next_index = 0;
    for (uint64_t i = 0; i < NumTargets; i++)
    {
        if (Dist[i] == 0) continue;
        sorted_dist[next_index++] = Dist[i];
    }
    sort(sorted_dist, sorted_dist + NumTargets - TargetsFound);

    // indices that can pass through the filter
    filter_dist = sorted_dist[min(nearest - 1, NumTargets-TargetsFound - 1)]; // largest distance that will pass through the filter
    for (uint64_t i = 0; i < NumTargets; i++)
    {

        if (Dist[i] <= filter_dist && Dist[i] > 0) 
        {
            filter_indices.push_back(i);
        }
    }

    // initialization
    uint64_t bestDist = LARGE;
    uint64_t currentDist;
    vector<uint64_t> candidates;


    for (uint64_t i = 0; i < counter; i++)
    {
        // Filtering
        if (!filtering(AllElements[i].newDist, filter_indices, Dist)) continue;
        // Normal BP rand
        currentDist = calculateDist(AllElements[i].newDist, NumTargets);
        if (currentDist < bestDist)
        {
            // Updating the best distance and norm
            bestDist = currentDist;
            // clear previous candidates
            candidates.clear();
            // inputting the new candidates
            candidates.push_back(i);
        }
        else if (currentDist == bestDist)
        {
            // equal candidates
            candidates.push_back(i);
        }
    }
    // randomly choose one of the candidates
    rand_generator.seed(time(0));
    uniform_int_distribution<int> rand_distribution(0, candidates.size() - 1);
    int rand_num = rand_distribution(rand_generator);
    return candidates[rand_num];
}

bool filtering(uint64_t tempDist[], vector<uint64_t> filter_indices, uint64_t *Dist)
{
    for (uint64_t i = 0; i < filter_indices.size(); i++)
    {
        // if any of the acceptable distance is reduced, return true immediately
        if (tempDist[filter_indices[i]] < Dist[filter_indices[i]]) return true;
    }
    return false;
}

uint64_t calculateDist(uint64_t A[],uint64_t length)
{
    uint64_t s = 0;
    for (uint64_t i = 0; i < length; i++) s += A[i];
    return s;
}

uint64_t calculateNorm(uint64_t A[],uint64_t length)
{
    uint64_t s = 0;
    for (uint64_t i = 0; i < length; i++) s += A[i] * A[i];
    return s;
}

void refreshDist(uint64_t NumTargets, uint64_t *Dist, uint64_t *InitDist)
{
    for(uint64_t k = 0; k < NumTargets; k++) Dist[k] = InitDist[k];
}


bool is_base(uint64_t x, uint64_t BaseSize, uint64_t *Base) 
{
    //sanity check, shouldn't ask if 0 is base
    if (x == 0) 
    {
        cout << "asking if 0 is in Base " << endl;
        exit(0);
    }

    for (uint64_t i = 0; i < BaseSize; i++) 
    {
        if (x == Base[i])
        {
            return true;
        }
    }
    return false;
}

uint64_t NewDistance(uint64_t u, uint64_t NumInputs, uint64_t BaseSize, uint64_t NewBase, uint64_t *Target, uint64_t *Base, uint64_t *Dist) 
{
	
    //if Target[u] is in augmented base return 0;
    if (is_base(Target[u], BaseSize, Base) || (NewBase == Target[u])) 
    {
        return 0;
    }

    if (reachable(Target[u] ^ NewBase, Dist[u] - 1, NumInputs, BaseSize, Base)) 
    {
        return (Dist[u] - 1);
    } 
    else
    {
        return Dist[u]; //keep old distance
    }
}

//return true if T is the sum of K elements among Base[S..BaseSize-1]
bool reachable(uint64_t T, uint64_t K, uint64_t S, uint64_t BaseSize, uint64_t *Base) 
{
    if (__builtin_popcount(T) <= K)
    {
        return true;
    }
    if (S > BaseSize - 1) 
    {
        return false; //exceeded count
    }

    if (K == 0) 
    {
        return false; //this is probably not reached
    }

    if (K == 1) 
    {
        for (int i = S; i < BaseSize; i++) if (T == Base[i]) 
        {
            return true;
        }
        return false;
    }
    
    //consider those sums containing Base[S]
    if (reachable(T ^ Base[S], K - 1, S + 1, BaseSize, Base)) 
    {
        return true;
    }

    //consider those sums not containing Base[S]
    if (reachable(T, K, S + 1, BaseSize, Base)) 
    {
        return true;
    }
    //not found
    return false;
}