#include"reduce.h"
#include"rnbp_a1_a2.h"
#include"paar.h"

void ReadTargetMatrix(uint64_t *NumInputs, uint64_t *NumTargets, uint64_t *Target, uint64_t *Dist, uint64_t *InitDist);
void read_implementation(vector<list> &seq);

int SIZE;
mt19937 rand_generator;

int main(int argc, char *argv[])
{	
	uint64_t NumInputs = 0;
    uint64_t NumTargets = 0;
	uint64_t Target[MaxBaseSize] = { 0 };
	uint64_t Dist[MaxBaseSize] = { 0 };
	uint64_t InitDist[MaxBaseSize] = { 0 };
	
	ReadTargetMatrix(&NumInputs, &NumTargets, Target, Dist, InitDist);
	vector<list> seq;
	
	if(FUNCTION == 0)	//FurOpt: Further optimise the given s-Xor implementation;
	{
		SIZE = NumInputs;
		read_implementation(seq);
		(NumTargets < NumInputs)? SIZE = NumInputs : SIZE = NumTargets;
	}

	if(FUNCTION == 1)	//IterOpt: Read the matrix, and search for its optimized implementation;
	{
 		rand_generator.seed(time(0));
    	uniform_int_distribution<int> rand_distribution(0, 5);
    	int OPTION_A = rand_distribution(rand_generator);
		(NumTargets < NumInputs)? SIZE = NumInputs : SIZE = NumTargets;
	
		for(int i = 0; i < SIZE; i++)
		{
			list tmp;
			tmp.value.push_back(i);
			tmp.flag = -1;
			tmp.usd = 0;
			seq.push_back(tmp); 
	 	}

    	switch(OPTION_A)
    	{
			case 0: Paar_algorithm1(seq, NumInputs, NumTargets, Target); break;
			case 1: RPaar_algorithm1(seq, NumInputs, NumTargets, Target); break;
			case 2: BP_algorithm(seq, NumInputs, NumTargets, Target, Dist); break;
			case 3: A1_A2_RNBP_algorithm(seq, NumInputs, NumTargets, 3, Target, Dist, InitDist); break;
			case 4: A1_A2_RNBP_algorithm(seq, NumInputs, NumTargets, 4, Target, Dist, InitDist); break;
			case 5: A1_A2_RNBP_algorithm(seq, NumInputs, NumTargets, 5, Target, Dist, InitDist); break;
			default:printf("error\n"); break;
    	}
	}
 	get_reduced(seq);

 	return 0;
}

void ReadTargetMatrix(uint64_t *NumInputs, uint64_t *NumTargets, uint64_t *Target, uint64_t *Dist, uint64_t *InitDist)
{	uint64_t N_in = 0;
	uint64_t N_out = 0;
    cin >> N_out;
    cin >> N_in;
    *NumInputs = N_in;
    *NumTargets = N_out;
    int bit;
    for (int i = 0; i < *NumTargets; i++) 
    {
        uint64_t PowerOfTwo  = 1;
        Target[i] = 0;
        Dist[i] = -1; //initial distance from Target[i] is Hamming weight - 1
        for (int j = 0; j < *NumInputs; j++) 
        {
            cin >> bit;
            if (bit) 
            {	
            	Dist[i]++; //after using randomly permutation, change the Dist[]!!!!!!;
                Target[i] = Target[i] + PowerOfTwo;
            }
            PowerOfTwo = PowerOfTwo << 1;
        }
    }
    for (int k = 0; k < *NumInputs; k++) InitDist[k] = Dist[k];
}

void write_file(vector<list> seq)
{
	ofstream f;
	f.open("Reduced_Result.txt");
	f << "Xor Count = " << seq.size() - SIZE << endl;
	for(int i = SIZE; i < seq.size(); i++)
	{ 
		f << "t["<<seq[i].value[0] << "] = t[" << seq[i].value[1] << "] ^ t[" << seq[i].value[2] << "]";
		if(seq[i].flag >= 0)
		{
			f << "    y[" << seq[i].flag << "]";
		}
		f << endl;
	}
	f.close();
}

//Read the given s-Xor implementation from Xiang et al..
void read_implementation(vector<list> &seq)
{
	fstream f("S-Xor-Implementation.txt");
	string line;

	for(int i = 0; i < SIZE; i++)
	{
		list tmp;
		tmp.value.push_back(i);
		tmp.flag = -1;
		tmp.usd = 0;
		seq.push_back(tmp); 
	 }

	while(getline(f, line))
	{
		list tmp_r;
		int num = 0;
		for(int i = 0; i < line.length(); i++)
		{
			if((line[i] >= '0') && (line[i] <= '9'))
			{
				num = num * 10 + (line[i] - '0');
			}
			else if(line[i] == ']')
			{
				tmp_r.value.push_back(num);
				num = 0;
			}
		}
		tmp_r.usd = 0;
		if(tmp_r.value.size() == 3)
		{
			tmp_r.flag = -1;
		}
		if(tmp_r.value.size() == 4)
		{
			tmp_r.flag = tmp_r.value[3];
			tmp_r.value.pop_back();
		}
		seq.push_back(tmp_r);
	}

	for(int i = SIZE; i < seq.size(); i++)
	{
		for(int j = i + 1; j < seq.size(); j++)
		{
			for(int k = 0; k < seq[i].value.size(); k++)
			{
				if(seq[j].value[k] == seq[i].value[0])
				{
					seq[j].value[k] = i;
				}
			}
		}
		seq[i].value[0] = i;
	}

	for(int i = SIZE; i < seq.size(); i++)
	{
		for(int j = 1; j< seq[i].value.size(); j++)
		{
			seq[seq[i].value[j]].usd += 1;
		}
	}
}
