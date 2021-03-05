#include <pthread.h>
#include <cstdlib>
#include"reduce.h"
#include"rnbp_a1_a2.h"

extern int SIZE;
int MAX_SIZE = 1000;

vector<uint64_t> Ori_Y;

void* reduce_thread(void *d)
{
	thread_data* data = (thread_data*)d;
	get_equivalent_circuit(data->seq, data->start, data->gap);
	pthread_exit(NULL);
}

#if MULTI_THREAD_FLAG
void reduce_gap(vector<list> &seq)
{
	int gap = seq.size() - SIZE - 1;//int gap = seq.size() - SIZE;
	thread_data data[THREAD_NUM];

	int s_p = SIZE;
	bool flag = true;
	while(flag)
	{
		for(int j = 0; j < THREAD_NUM; j++)
		{
			data[j].seq = seq;
			data[j].gap = gap;
			data[j].start = s_p;
			//cout << "gap = " << gap << endl;
			if((gap != 1)||(s_p != (seq.size() - 1)))
			{
				if(s_p == (seq.size() - gap))
				{
					gap--;
					s_p = SIZE;
				}
				else
				{
					s_p++;
				}
			}
			else
			{
				flag = false;
			}
		}

		int rc;
		pthread_t threads[THREAD_NUM];
		pthread_attr_t attr;
		void* status;
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
		for(int j = 0; j < THREAD_NUM; j++)
		{
			rc = pthread_create(&threads[j], NULL, reduce_thread, (void *)(data + j));
			if(rc)
			{
				cout << "Error: unable to create thread!" << endl;
				exit(-1);
			}
		}
		// pthread_attr_destroy(&attr);
		for(int j = 0; j < THREAD_NUM; j++)
		{
			rc = pthread_join(threads[j], &status);
			if(rc)
			{
				cout << "Error: unable to join, " << rc << endl;
				exit(-1);
			}
		}
		for(int j = 0; j < THREAD_NUM; j++)
		{
			while(findreduce(data[j].seq))
			{
				update_seq(data[j].seq);
			}

			if(data[j].seq.size() < seq.size())
			{
				seq = data[j].seq;
				gap = data[j].seq.size() - SIZE - 1;
				s_p = SIZE;
				if((check(seq)) && ((seq.size() - SIZE) < MAX_SIZE))
	    		{	
	    			MAX_SIZE = seq.size() - SIZE;
					write_file(seq);
					cout << "-------------------------Optimal seq find at present: seq.size() = " << seq.size() - SIZE << endl;
				}
			}
		}
	}
}

#else
void reduce_gap(vector<list> &seq)
{
	int gap = seq.size() - SIZE - 1;
	while(gap >= 2)
	{   
		int bf_len = seq.size();
		int s_p = SIZE;
        	
		while((s_p + gap) <= seq.size())
		{  
			cout << "gap = " << gap << ";    " << "start_p = " << s_p << ";    -----Optimised seq.size =  " << MAX_SIZE << endl;
			vector<list> t_seq;
			t_seq.assign(seq.begin(), seq.end());
			get_equivalent_circuit(seq, s_p, gap);
				
			while(findreduce(seq))
			{
				update_seq(seq);
			}

			if(seq.size() < t_seq.size())
			{  
				if((check(seq)) && ((seq.size() - SIZE) < MAX_SIZE))
	    		{	
	    			MAX_SIZE = seq.size() - SIZE;
					write_file(seq);
					cout << "-------------------------Optimal seq find at present: seq.size() = " << MAX_SIZE << endl;
				}
				gap = seq.size() - SIZE - 1;//seq.size() - SIZE
				s_p = SIZE;
			}
        	else 
			{
				seq.assign(t_seq.begin(), t_seq.end());
				s_p++;
			}
		}	
        	
		if(seq.size() < bf_len)
        {
        	gap = seq.size() - SIZE - 1;
		}
		else if(seq.size() >= bf_len)
		{
			gap--;
		}
    }
}

#endif
int get_reduced(vector<list> &o_seq)
{
	cout << "Original seq.size() = " << o_seq.size() - SIZE << endl;
	Ori_Y = find_targets(o_seq);
	update_seq(o_seq);
	while(findreduce(o_seq))
	{  
		update_seq(o_seq);
	}
	if((check(o_seq)) && ((o_seq.size() - SIZE) < MAX_SIZE))
	{	
	    MAX_SIZE = o_seq.size() - SIZE;
		write_file(o_seq);
		cout << "Original seq.size() = " << o_seq.size() - SIZE << ";     XorCount found  = " << MAX_SIZE << endl;
	}
	/////////////////////
	while(1)
	{  
		vector<list> seq(o_seq);
		if(FUNCTION == 1)
		{
			get_equivalent_circuit(seq, SIZE, o_seq.size() - SIZE);
			cout << "Original seq.size() = " << seq.size() - SIZE << ";     XorCount found  = " << MAX_SIZE << endl;
			while(findreduce(seq))
			{  
				update_seq(seq);
			}
		
	   		if((check(seq)) && ((seq.size() - SIZE) < MAX_SIZE))
	    	{	
	    		MAX_SIZE = seq.size() - SIZE;
				write_file(seq);
				cout << "----------------------- " << "    XorCount found  = " << MAX_SIZE << endl;
			}
		}
		reduce_gap(seq);
	}

	return 0;
}

bool check(vector<list> seq)
{
	vector<uint64_t> y = find_targets(seq);
	for(int i = 0; i < SIZE; i++)
	{
		if(y[i] != Ori_Y[i])
		{
			cout << "targets wrong at:  " << i << endl;
			return false;
		}
	}

	return true;
}

///////get equivalent seq. from algorithm B////////
void get_equivalent_circuit(vector<list> &seq, int sp, int length)
{
	vector<list> tmp_seq;
	vector<int> input;
	//expand the segment chosen
	pre_bp(seq, sp, length);

	for(int i = SIZE; i < seq.size(); i++)
	{   
		delete_twice_ele(seq[i].value);
		sort(seq[i].value.begin() + 1, seq[i].value.end());
	}

	while(pre_reduce_seq(seq))
	{
		for(int i = SIZE; i < seq.size(); i++)
		{
			sort(seq[i].value.begin() + 1, seq[i].value.end());
		}
	}
	
	int dst = 0;
	for(int i = sp; i < seq.size(); i++)
	{
		if(seq[i].value[0] < length + sp)
		{
			dst++;
		}
	}

	count_used_time(seq);

	tmp_seq.assign(seq.begin() + sp, seq.begin() + sp + dst); 
	//the length of chosen segment has changed from 'gap' to 'dst';
	//delete the register never been used;
	for(int i = tmp_seq.size() - 1; i >= 0; i--)
	{   
		if((tmp_seq[i].usd == 0) && (tmp_seq[i].flag < 0))
		{
			tmp_seq.erase(tmp_seq.begin() + i);
		}
	}
	
	//find the real input/output bits;
	uint64_t Target[MaxBaseSize] = { 0 };
	uint64_t Dist[MaxBaseSize] = { 0 };
	uint64_t InitDist[MaxBaseSize] = { 0 };
	find_inputs(tmp_seq, input);
	find_outputs(tmp_seq, input, Target, Dist, InitDist);

	vector<list> e_seq = Algorithm_B(tmp_seq, input, Target, Dist, InitDist);
	seq.erase(seq.begin() + sp, seq.begin() + sp + dst);

	int sp_B = input.size();
	seq.insert(seq.begin() + sp, e_seq.begin() + sp_B, e_seq.end());
	update_seq(seq);
}

void find_inputs(vector<list> tmp_seq, vector<int> &input)
{
	for(int i = 0; i < tmp_seq.size(); i++)
	{	
		input.insert(input.end(), tmp_seq[i].value.begin() + 1, tmp_seq[i].value.end());
	}
	sort(input.begin(), input.end());
	input.erase(unique(input.begin(), input.end()), input.end());
}

void find_outputs(vector<list> tmp_seq, vector<int> input, uint64_t *Target, uint64_t *Dist, uint64_t *InitDist)
{   	
	uint64_t NumTargets = tmp_seq.size();
    uint64_t NumInputs = input.size();
    int SIZE_B = NumTargets;
    if(NumTargets < NumInputs)
    {
        SIZE_B = NumInputs;
    }
    /*//check that NumInputs is < wordsize
    if (NumInputs >= 8*sizeof(uint64_t)) 
    {
        cout << "        too many inputs" << endl;
		return false;
        //exit(0);
    }*/

    for (uint64_t i = 0; i < NumTargets; i++) 
    { //read row i
        Target[i] = 0;
        Dist[i] = -1; //initial distance from Target[i] is Hamming weight - 1
        for (uint64_t k = 1; k < tmp_seq[i].value.size(); k++)
        {
            uint64_t PowerOfTwo  = 1;
            for(uint64_t j = 0; j < input.size(); j++)
            {
                if(input[j] == tmp_seq[i].value[k])
                {
                    Dist[i]++;
                    Target[i] = Target[i] + PowerOfTwo;
                }
                PowerOfTwo = PowerOfTwo << 1;
            }    
        }
    }
    // Update the InitDist for subsequent rounds
	//cout<< "        Matrix into BP is : " << tmp_seq.size() << " x " << input.size() <<endl;
    for (uint64_t k = 0; k < NumTargets; k++) InitDist[k] = Dist[k];
}

//expand the rest registers
void pre_bp(vector<list> &seq, int sp, int gap)
{
	vector<list> tmp_seq;
	tmp_seq.insert(tmp_seq.begin(), seq.begin(), seq.end());
	for(int i = sp; i < sp + gap; i++)
	{   
		delete_twice_ele(tmp_seq[i].value);
		sort(tmp_seq[i].value.begin() + 1, tmp_seq[i].value.end());
	}

	for(int i = SIZE; i < sp + gap - 1; i++)
	{
		for(int j = i + 1; j < sp + gap; j++)
		{
			for(int k = 1; k < tmp_seq[j].value.size(); k++)
			{
				if(tmp_seq[j].value[k] == tmp_seq[i].value[0])
				{
					tmp_seq[j].value.erase(tmp_seq[j].value.begin() + k);
					for(int p = 1; p < tmp_seq[i].value.size(); p++)
					{   
						vector<int>::iterator iter = find(tmp_seq[j].value.begin() + 1, tmp_seq[j].value.end(), tmp_seq[i].value[p]);
						if(iter != tmp_seq[j].value.end())
						{
							iter = tmp_seq[j].value.erase(iter);
						}
						else
						{
							tmp_seq[j].value.push_back(tmp_seq[i].value[p]);
						}
					}  
				} 
			}	
		}
	}
	seq.erase(seq.begin() + sp, seq.begin() + sp + gap);
	seq.insert(seq.begin() + sp, tmp_seq.begin() + sp, tmp_seq.begin() + sp + gap);
}

//for the case that : t_i = t_j ^ t_j; i.e., t_i = NULL, storing nothing and has never been used in seq;
bool delete_r_useless(vector<list> &seq)
{
	for(int i =  SIZE; i < seq.size(); i++)
	{
		if((seq[i].usd == 0) && (seq[i].flag < 0))
		{
			seq.erase(seq.begin() + i);
			return true;
		}
	}
	return false;
}

//for the case that : t_i = t_j ^ t_j; i.e., t_i = NULL, storing nothing but has been used in seq;
bool delete_r_null(vector<list> &seq, int idx)
{
	if((seq[idx].value.size() == 1) && (seq[idx].flag < 0))//if(seq[i].value.size() == 1)
	{   
		for(int j = SIZE; j < seq.size(); j++)
		{
			for(int k = seq[j].value.size() - 1; k >= 1 ; k--)
			{
				if(seq[j].value[k] == seq[idx].value[0])
				{
					seq[j].value.erase(seq[j].value.begin() + k);
				}
			}
		}
		seq.erase(seq.begin() + idx);
		return true;
	}
	return false;
}

bool delete_r_replace(vector<list> &seq, int idx)
{
	if((seq[idx].value.size() == 2) && (seq[idx].value[0] != seq[idx].value[1]))
	{   
		for(int j = SIZE; j < seq.size(); j++)
		{	//replace t10 by t2 in seq;
			for(int k = 1; k < seq[j].value.size(); k++)
			{
				if(seq[j].value[k] == seq[idx].value[0])
				{
					seq[j].value[k] = seq[idx].value[1];
				}
			}
		}
		//if t10 = y; let t2 = y; and then, delete t10 in seq; 
		if(seq[idx].flag >= 0)
		{
			for(int k = SIZE; k < seq.size(); k++)
			{
				if(seq[k].value[0] == seq[idx].value[1])
				{
					seq[k].flag = seq[idx].flag;
				}
			}
		}
		seq.erase(seq.begin() + idx);
		return true;
	}
	return false;
}

bool delete_r_same(vector<list> &seq, int idx)
{
	for(int j = seq.size() - 1; j > idx; j--)
	{
		if(seq[idx].value.size() == seq[j].value.size())
		{
			bool equal = true;
			for(int k = 1; k < seq[idx].value.size(); k++)
			{
				if(seq[idx].value[k] != seq[j].value[k])
				{
					equal = false;
				}
			}
			if(equal)
			{	//t10 = t0 + t2（ + t3 + ...); t20 = t0 + t2（ + t3 + ...); delete t20;
				for(int k = SIZE; k < seq.size(); k++)
				{	//replace t20 by t10 in seq;
					for(int p = 1; p < seq[k].value.size(); p++)
					{
						if(seq[k].value[p] == seq[j].value[0])
						{
							seq[k].value[p] = seq[idx].value[0];
						}
					}
				}
				if(seq[j].flag >= 0)
			  	{
					seq[idx].flag = seq[j].flag;	
				}
				seq.erase(seq.begin() + j);
				return true;
			}
		}
	}
	return false;
}

bool pre_reduce_seq(vector<list> &seq)
{  
	if(delete_r_useless(seq))
	{
		return true;
	}
	for(int i =  SIZE; i < seq.size(); i++)
	{	// t10 = t0 + t10 ( + ti + ...); if t_i include itself, delete t_i in it;
		vector<int>::iterator iter = find(seq[i].value.begin() + 1, seq[i].value.end(), seq[i].value[0]);
		if(iter != seq[i].value.end())
		{
			iter = seq[i].value.erase(iter);
			return true;
		}
		else
		{   
			// t10 = null; delete t10 in each register t_i that include t10; delete t10 in seq;
			if(delete_r_null(seq, i))
			{
				return true;
			}
			// t10 = t2;
			if(delete_r_replace(seq, i))
			{
				return true;
			}
			//if two registers storing the same elements: delete the later one;
			if(delete_r_same(seq, i))
			{
				return true;
			}
		}
	}
	return false;
}

bool delete_twice_ele(vector<int> &aray)
{
	for(int i = 1; i < aray.size(); i++)
    {
      for(int j = aray.size() - 1; j > i; j--)
      {
        if(aray[i] == aray[j])
        {
          aray.erase(aray.begin() + j);
          aray.erase(aray.begin() + i);
          return true;
        }
      }
    }
    return false;
}