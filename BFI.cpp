#include"bfi.h"
vector<uint64_t> generate_p(int number)
{
	vector<uint64_t> B;
	vector<uint64_t> P;
	for(int j = 0; j < number; j++)
	{
		B.push_back(j);
	}
	
	while(B.size() > 0)
	{	
		uint64_t a = rand() % B.size();
		vector<uint64_t>::iterator iter = find(B.begin(), B.end(), B[a]);
		if(iter != B.end())
		{	
			P.push_back(B[a]);
			iter = B.erase(iter);
		}
		else
		{
			cout << "generate random permutation wrong!" << endl;
		}
	}
	
	return P;
}

vector<uint64_t> p_m_q(uint64_t *Target, uint64_t no_row, uint64_t no_col, uint64_t *Dist, vector<uint64_t> P, vector<uint64_t> Q)
{		
	vector<uint64_t> target_bar;
	vector<uint64_t> target_reuslt;

	uint64_t tmp_dis[no_row] = { 0 };
	//update Dist
	for(int i = 0; i < no_row; i++)
	{
		tmp_dis[i] = Dist[P[i]];
	}
	for(int i = 0; i < no_row; i++)
	{
		Dist[i] = tmp_dis[i];
	}

	for(int i = 0; i < no_row; i++)
	{	
		uint64_t sum = 0;
		for(uint64_t j = 0; j < no_col; j++)
		{	
			sum = sum + (((Target[i] >> Q[j]) & 1ull) << j);
		}
		target_bar.push_back(sum);
	}

	for(int j = 0; j < no_row; j++)
	{
		target_reuslt.push_back(target_bar[P[j]]);
	}

	for(int j = 0; j < no_row; j++)
	{
		Target[j] = target_reuslt[j];
	}

	return target_reuslt;
}

vector<uint64_t> p_m_q_paar(uint64_t *Target, uint64_t no_row, uint64_t no_col, vector<uint64_t> P, vector<uint64_t> Q)
{		
	vector<uint64_t> target_bar;
	vector<uint64_t> target_reuslt;

	for(int i = 0; i < no_row; i++)
	{	
		uint64_t sum = 0;
		for(uint64_t j = 0; j < no_col; j++)
		{	
			sum = sum + (((Target[i] >> Q[j]) & 1ull) << j);
		}
		target_bar.push_back(sum);
	}

	for(int j = 0; j < no_row; j++)
	{
		target_reuslt.push_back(target_bar[P[j]]);
	}

	for(int j = 0; j < no_row; j++)
	{
		Target[j] = target_reuslt[j];
	}

	return target_reuslt;
}

void recover_implementation(vector<list> &seq, vector<uint64_t> P, vector<uint64_t> Q)
{
	vector<list> tmp_seq;
	tmp_seq.assign(seq.begin(), seq.end());
	
	for(int i = Q.size(); i < seq.size(); i++)
	{	
		for(int j = 1; j < seq[i].value.size(); j++)
		{
			if(seq[i].value[j] < Q.size())
			{
				seq[i].value[j] = Q[tmp_seq[i].value[j]];
			}
		}
		if(seq[i].flag >= 0)
		{
			seq[i].flag = P[tmp_seq[i].flag];
		}
	}
}
