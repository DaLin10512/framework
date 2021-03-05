#include"paar.h"
#include"rnbp_a1_a2.h"

extern mt19937 rand_generator;
extern int SIZE; 

int hamming_weight(uint64_t input) 
{
    return __builtin_popcount(input & 0xffffffff) + __builtin_popcount((input >> 32) & 0xffffffff);
}

int Paar_algorithm1(vector<list> &seq, uint64_t NumInputs, uint64_t NumTargets, uint64_t *Target) 
{   
    uint64_t *input_matrix;
    input_matrix = (uint64_t *)malloc((SIZE + 1000) * sizeof(uint64_t));

    vector<uint64_t> P;
    vector<uint64_t> Q;
   	int randomly_flag = rand() % 2;
   	if(randomly_flag == 1)
    {
    	P = generate_p(NumTargets);
    	Q = generate_p(NumInputs);
    	p_m_q_paar(Target, NumTargets, NumInputs, P, Q);
    }

    for(int i = 0; i < NumInputs; i++)
    {  
        input_matrix[i] = 0;
        for(int j = 0; j < NumTargets; j++)
        {
            input_matrix[i] += ((Target[j] >> i) & 1ull) << j;
        }
    }

    int xor_count = 0; /*In fact, here we do not need to record the xor_count, since 'seq.size() - SIZE' is equal to it.*/
    uint64_t number_of_columns = NumInputs;
    int hw_max = 0;
    int i_max = 0, j_max = 0;
    uint64_t tmp = 0;
    int hw = 0;
    uint64_t new_column = 0;
    vector<pair<int, int>> program;

    do {
        hw_max = 0;
        for (uint64_t i = 0; i < number_of_columns; i++) 
        {  
            for (uint64_t j = i + 1; j < number_of_columns; j++) 
            {  
                tmp = input_matrix[i] & input_matrix[j];
                hw = hamming_weight(tmp);
                if (hw > hw_max) 
                {
                    hw_max = hw;
                    i_max = i;
                    j_max = j;
                }
            }
        }
        if (hw_max >= 1) 
        {
            new_column = input_matrix[i_max] & input_matrix[j_max];
            input_matrix[number_of_columns] = new_column;
            input_matrix[i_max] =  new_column ^ input_matrix[i_max];
            input_matrix[j_max] =  new_column ^ input_matrix[j_max];
            xor_count -= (hw_max - 1);
            number_of_columns++;
            program.push_back(make_pair(i_max, j_max));
        }
    } while (hw_max >= 1);

    int ctr = NumInputs;
    
    for(auto const& prog:program) 
    {   
        list row;
        row.flag = -1;
        row.usd = 0;
        row.value.push_back(ctr);
        row.value.push_back(prog.first);
        row.value.push_back(prog.second);
        seq.push_back(row);
        ctr++;
    }
  
    for(uint64_t i = 0; i < NumTargets; i++)
    {
        for(uint64_t j = 0; j < number_of_columns; j++) 
        {
            if ((input_matrix[j] & (1ll << (NumTargets - 1 - i))) != 0) 
            {
                for(int k = 0; k < seq.size(); k++)
                {
                    if(seq[k].value[0] == j)
                    {
                        seq[k].flag = (NumTargets - 1 - i);
                    }
                }
            }
        }
    }
    free(input_matrix);

    if(randomly_flag == 1)
    {
    	recover_implementation(seq, P, Q);
	}
 
    return xor_count;
}

int RPaar_algorithm1(vector<list> &seq, uint64_t NumInputs, uint64_t NumTargets, uint64_t *Target) 
{   
	uint64_t *input_matrix;
    input_matrix = (uint64_t *)malloc((SIZE + 1000) * sizeof(uint64_t));

    vector<uint64_t> P;
    vector<uint64_t> Q;
    int randomly_flag = rand() % 2;
	if(randomly_flag == 1)
    {
    	P = generate_p(NumTargets);
    	Q = generate_p(NumInputs);
    	p_m_q_paar(Target, NumTargets, NumInputs, P, Q);
	}

    for(int i = 0; i < NumInputs; i++)
    {  
        input_matrix[i] = 0;
        for(int j = 0; j < NumTargets; j++)
        {
            input_matrix[i] += ((Target[j] >> i) & 1) << j;
        }
    }

    int xor_count = 0;
    uint64_t number_of_columns = NumInputs;
    int hw_max = 0;;
    int i_max = 0, j_max = 0;
    uint64_t tmp = 0;
    uint64_t new_column = 0;
    // compute naive xor count
    for (int i = 0; i < NumInputs; i++) 
    {
        xor_count += hamming_weight(input_matrix[i]);
    }
    xor_count -= NumTargets;
    vector<pair<int, int>> program;

    do {
        vector<pair<int, int>> idx_i_j = find_max_idx(input_matrix, number_of_columns);
        i_max = idx_i_j[0].first; 
        j_max = idx_i_j[0].second;
        tmp = input_matrix[i_max] & input_matrix[j_max];
        hw_max = hamming_weight(tmp);
        if (hw_max >= 1) 
        {
            new_column = input_matrix[i_max] & input_matrix[j_max];
            input_matrix[number_of_columns] = new_column;
            input_matrix[i_max] = new_column ^ input_matrix[i_max];
            input_matrix[j_max] = new_column ^ input_matrix[j_max];
            xor_count -= (hw_max - 1);
            number_of_columns++;
            program.push_back(make_pair(i_max, j_max));
        }
    } while (hw_max >= 1);

    int ctr = NumInputs;

    for(auto const& prog:program) 
    {   
        list row;
        row.flag = -1;
        row.usd = 0;
        row.value.push_back(ctr);
        row.value.push_back(prog.first);
        row.value.push_back(prog.second);
        seq.push_back(row);
        //cout << "x" << ctr << " = x" << prog.first << " + x" << prog.second << endl;
        ctr++;
    }
   
    for(uint64_t i = 0; i < NumTargets; i++)
    {
        for(uint64_t j = 0; j < number_of_columns; j++) 
        {
            if ((input_matrix[j] & (1ll << (NumTargets - 1 - i))) != 0) 
            {
                //cout << " = x" << j;
                for(int k = 0; k < seq.size(); k++)
                {
                    if(seq[k].value[0] == j)
                    {
                        seq[k].flag = (NumTargets - 1 - i);
                    }
                }
            }
        }
    }
    free(input_matrix);

    if(randomly_flag == 1)
    {
    	recover_implementation(seq, P, Q);
	}
    return xor_count;
}

vector<pair<int, int>> find_max_idx(uint64_t *input_matrix, int number_of_columns)
{
    uint64_t tmp;
    int hw = 0;
    int i_max_list[10000];
    int j_max_list[10000];
    int hw_max = 0;
    int list_end = 0; 
    for (int i = 0; i < number_of_columns; i++) 
    {
        for (int j = i + 1; j < number_of_columns; j++) 
        {
            tmp = input_matrix[i] & input_matrix[j];
            hw = hamming_weight(tmp);
            if (hw > hw_max) {
                hw_max = hw;
                i_max_list[0] = i;
                j_max_list[0] = j;
                list_end = 1;
            }
       
            else if ((hw == hw_max)&& (hw != 0)) 
            {
                i_max_list[list_end] = i;
                j_max_list[list_end] = j;
                list_end++;
            }
        }
    }
    
    vector<pair<int, int>> idx_i_j_max;
    int rand_num = 0;
    if(list_end > 1)
    {
        rand_generator.seed(time(0));
        uniform_int_distribution<int> rand_distribution(0, list_end - 1);
        rand_num = rand_distribution(rand_generator);
    }
    idx_i_j_max.push_back(make_pair(i_max_list[rand_num], j_max_list[rand_num]));

    return idx_i_j_max;
}
