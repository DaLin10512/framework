#include"paar.h"
#include"regidster.h"
extern int SIZE;
extern int Check_NumInputs;
extern int Check_NumTargets;

vector<uint64_t> find_targets(vector<list> seq)
{
  vector<uint64_t> check_y;
	vector<uint64_t> c_tmp_y;
	
	for(int i = 0; i < SIZE; i++)
	{
    	c_tmp_y.push_back(1ull << i);//for linux   
	}
	for(int i = SIZE; i < seq.size(); i++)
	{
		uint64_t tmp = 0;
		tmp = c_tmp_y[seq[i].value[1]] ^ c_tmp_y[seq[i].value[2]];
		c_tmp_y.push_back(tmp);
	}
	int idx = 0;
	while(check_y.size() < SIZE)
	{
		for(int i = SIZE; i < seq.size(); i++)
		{
			if(seq[i].flag == idx)
			{
				check_y.push_back(c_tmp_y[i]);
				idx++;
			}
		}
	}

	return check_y;
}

void update_seq(vector<list> &seq)
{
	sort_seq(seq);
	count_used_time(seq);
	while(pre_reduce_seq(seq))
	{
		sort_seq(seq);
		count_used_time(seq);
	}	
}

void sort_seq(vector<list> &seq)
{   
	for(int i = SIZE; i < seq.size(); i++)
	{
		if(seq[i].value[0] != i + 100000)
		{   
			for(int j = SIZE; j < seq.size(); j++) 
			{
				for(int k = 1; k < seq[j].value.size(); k++)
				{  
					if(seq[j].value[k] == seq[i].value[0])
					{
						seq[j].value[k] = i + 100000;
					}
				}
			}
			seq[i].value[0] = i + 100000;
		}    
	}

	for(int i = SIZE; i < seq.size(); i++)
	{
		for(int j = 0; j < seq[i].value.size(); j++)
		{
			if(seq[i].value[j] > 100000)
			{
				seq[i].value[j] = seq[i].value[j] - 100000;
			}
		}
	}
}

void count_used_time(vector<list> &seq)
{
  for(int i = SIZE; i < seq.size(); i++)
	{
    int num = 0;
		for(int j = i + 1; j < seq.size(); j++)
		{
			for(int k = 1; k < seq[j].value.size(); k++)
			{
        		if(seq[j].value[k] == seq[i].value[0]){num ++;}
			}
		}
		seq[i].usd = num;
	}
}

bool redece_rule_1(vector<list> &seq, int u, int v)
{
  if((seq[v].value[1] == seq[u].value[0]) && (seq[v].value[2] == seq[u].value[1]))
  {
      rule_1(seq, u, v, 2);
      return true;
  }

  if((seq[v].value[1] == seq[u].value[0]) && (seq[v].value[2] == seq[u].value[2]))
  {
      rule_1(seq, u, v, 1);
      return true;
  } 
    
  if((seq[v].value[2] == seq[u].value[0]) && (seq[v].value[1] == seq[u].value[1]))
  {
      rule_1(seq, u, v, 2);
      return true;
  } 
      
  if((seq[v].value[2] == seq[u].value[0]) && (seq[v].value[1] == seq[u].value[2]))
  {
      rule_1(seq, u, v, 1);
      return true;
  }
  return false;       
}

void rule_1(vector<list> &seq, int u, int v, int ind_u)
{
  for(int i = SIZE; i < seq.size(); i++)//start from SIZE
  {
    for(int j = 1; j < seq[i].value.size(); j++)
    { 
      if(seq[i].value[j] == seq[v].value[0])
      {
        seq[i].value[j] = seq[u].value[ind_u];
      }
    }
  }
               
  if(seq[v].flag >= 0)
  {
    seq[seq[u].value[ind_u]].flag = seq[v].flag;  
  }
  seq.erase(seq.begin() + v);//delete_ele(seq, v);
  //if u is only used for v = u + a; u < v, so, delete v, and then, delete u;
  if((seq[u].flag < 0) && ((seq[u].usd <= 1)))
  {
    seq.erase(seq.begin() + u);//delete_ele(seq, u);
  }
} 
    
//t[u]==a+b, t[v]==c+u, t[w]==a(b)+c   
bool reduce_rule_2(vector<list> &seq, int u, int v, int w)
{ 
  if((seq[v].value[1] == seq[u].value[0]) && (seq[v].value[2] != seq[u].value[1]) && (seq[v].value[2] != seq[u].value[2]))
  {
    if((seq[u].usd <= 1) && (seq[u].flag < 0))
    {
      if((seq[v].value[2] == seq[w].value[2]) && (seq[w].value[1] == seq[u].value[1]))
      {
        rule_2(seq, u, v, w, 2);
        return true;
      }
      if((seq[v].value[2] == seq[w].value[2]) && (seq[w].value[1] == seq[u].value[2]))
      {
        rule_2(seq, u, v, w, 1);
        return true;
      }
      if((seq[v].value[2] == seq[w].value[1]) && (seq[w].value[2] == seq[u].value[1]))
      {
        rule_2(seq, u, v, w, 2);
        return true;
      }
      if((seq[v].value[2] == seq[w].value[1]) && (seq[w].value[2] == seq[u].value[2]))
      {
        rule_2(seq, u, v, w, 1);
        return true;
      }
    }
  }

  if((seq[v].value[2] == seq[u].value[0]) && (seq[v].value[1] != seq[u].value[1]) && (seq[v].value[1] != seq[u].value[2]))
  {
    if((seq[u].usd <= 1) && (seq[u].flag < 0))
    {
      if((seq[v].value[1] == seq[w].value[1]) && (seq[w].value[2] == seq[u].value[1]))
      {
        rule_2(seq, u, v, w, 2);
        return true;
      }
      if((seq[v].value[1] == seq[w].value[1]) && (seq[w].value[2] == seq[u].value[2]))
      {
        rule_2(seq, u, v, w, 1);
        return true;
      }
      if((seq[v].value[1] == seq[w].value[2]) && (seq[w].value[1] == seq[u].value[1]))
      {
        rule_2(seq, u, v, w, 2);
        return true;
      }
      if((seq[v].value[1] == seq[w].value[2]) && (seq[w].value[1] == seq[u].value[2]))
      {
        rule_2(seq, u, v, w, 1);
        return true;
      }
    }
  }

  return false;
}

void rule_2(vector<list> &seq, int u, int v, int w, int idx_u)
{ 
  seq[v].value[2] = seq[u].value[idx_u];
  seq[v].value[1] = seq[w].value[0];
  seq.insert(seq.begin() + v, seq[w]);//insert_list(seq, v, w);
  seq.erase(seq.begin() + w + 1);//delete_ele(seq, w + 1); 
  seq.erase(seq.begin() + u);//delete_ele(seq, u);
}

//t[u]==a+b, t[v]==c+b, t[w]==u+f 
bool reduce_rule_6(vector<list> &seq, int u, int v, int w)
{
  if((seq[w].value[1] == seq[u].value[0]) && (seq[w].value[2] != seq[u].value[1]) && (seq[w].value[2] != seq[u].value[2]))
  {
    if((seq[u].usd <= 1) && (seq[u].flag < 0))
    {
      if((seq[w].value[2] == seq[v].value[2]) && (seq[u].value[2] == seq[v].value[1]))
      {
        rule_6(seq, u, v, w, 1);
        return true;
      }
      if((seq[w].value[2] == seq[v].value[1]) && (seq[u].value[2] == seq[v].value[2]))
      {
        rule_6(seq, u, v, w, 1);
        return true;
      }
      if((seq[w].value[2] == seq[v].value[2]) && (seq[u].value[1] == seq[v].value[1]))
      {
        rule_6(seq, u, v, w, 2);
        return true;
      }
      if((seq[w].value[2] == seq[v].value[1]) && (seq[u].value[1] == seq[v].value[2]))
      {
        rule_6(seq, u, v, w, 2);
        return true;
      }
    }
  } 

  if((seq[w].value[2] == seq[u].value[0]) && (seq[w].value[1] != seq[u].value[1]) && (seq[w].value[1] != seq[u].value[2]))
  {
    if((seq[u].usd <= 1) && (seq[u].flag < 0))
    {
      if((seq[w].value[1] == seq[v].value[2]) && (seq[u].value[2] == seq[v].value[1]))
      {
        rule_6(seq, u, v, w, 1);
        return true;
      }
      if((seq[w].value[1] == seq[v].value[1]) && (seq[u].value[2] == seq[v].value[2]))
      {
        rule_6(seq, u, v, w, 1);
        return true;
      }
      if((seq[w].value[1] == seq[v].value[2]) && (seq[u].value[1] == seq[v].value[1]))
      {
        rule_6(seq, u, v, w, 2);
        return true;
      }
      if((seq[w].value[1] == seq[v].value[1]) && (seq[u].value[1] == seq[v].value[2]))
      { 
        rule_6(seq, u, v, w, 2);
        return true;
      }
    }
  }

  return false;
}
    
void rule_6(vector<list> &seq, int u, int v, int w, int idx_u)
{ 
  seq[w].value[2] = seq[u].value[idx_u];
  seq[w].value[1] = seq[v].value[0];
  seq.erase(seq.begin() + u);//delete_ele(seq, u);
}

//t[u]==a+b, t[v]==c+d, t[w]==v+f   
bool reduce_rule_3(vector<list> &seq, int u, int v, int w)
{
  if((seq[w].value[1] == seq[v].value[0]) && (seq[w].value[2] != seq[v].value[1]) && (seq[w].value[2] != seq[v].value[2]))
  {
    if((seq[v].usd <= 1) && (seq[v].flag < 0))
    {
      if((seq[w].value[2] == seq[u].value[1]) && (seq[u].value[2] == seq[v].value[1]))
      {
        rule_3(seq, u, v, w, 2);
        return true;
      }
      if((seq[w].value[2] == seq[u].value[1]) && (seq[u].value[2] == seq[v].value[2]))
      {
        rule_3(seq, u, v, w, 1);
        return true;
      }
      if((seq[w].value[2] == seq[u].value[2]) && (seq[u].value[1] == seq[v].value[1]))
      {
        rule_3(seq, u, v, w, 2);
        return true;
      }
      if((seq[w].value[2] == seq[u].value[2]) && (seq[u].value[1] == seq[v].value[2]))
      {
        rule_3(seq, u, v, w, 1);
        return true;
      }
    }
  } 

  if((seq[w].value[2] == seq[v].value[0]) && (seq[w].value[1] != seq[v].value[1]) && (seq[w].value[1] != seq[v].value[2]))
  {
    if((seq[v].usd <= 1) && (seq[v].flag < 0))
    {
      if((seq[w].value[1] == seq[u].value[1]) && (seq[u].value[2] == seq[v].value[1]))
      {
        rule_3(seq, u, v, w, 2);
        return true;
      }
      if((seq[w].value[1] == seq[u].value[1]) && (seq[u].value[2] == seq[v].value[2]))
      {
        rule_3(seq, u, v, w, 1);
        return true;
      }
      if((seq[w].value[1] == seq[u].value[2]) && (seq[u].value[1] == seq[v].value[1]))
      {
        rule_3(seq, u, v, w, 2);
        return true;
      }
      if((seq[w].value[1] == seq[u].value[2]) && (seq[u].value[1] == seq[v].value[2]))
      { 
        rule_3(seq, u, v, w, 1);
        return true;
      }
    }
  }

  return false;
}
    
void rule_3(vector<list> &seq, int u, int v, int w, int idx_v)
{ 
  seq[w].value[2] = seq[v].value[idx_v];
  seq[w].value[1] = seq[u].value[0];
  seq.erase(seq.begin() + v);//delete_ele(seq, v);
}

//t[u]==a+b, t[v]==u+d, t[w]==v+f   
bool reduce_rule_4(vector<list> &seq, int u, int v, int w)
{
  if((seq[v].value[1] == seq[u].value[0]) && (seq[v].value[2] != seq[u].value[1]) && (seq[v].value[2] != seq[u].value[2]))
  {
    if((seq[u].usd <= 1) || (seq[v].usd <= 1))
    {
      if((seq[w].value[1] == seq[v].value[0]) && (seq[w].value[2] == seq[u].value[1]))
      {
        if(rule_4(seq, u, v, w, 2, 2))
        {
          return true;
        }
      }
      if((seq[w].value[1] == seq[v].value[0]) && (seq[w].value[2] == seq[u].value[2]))
      {
        if(rule_4(seq, u, v, w, 1, 2))
        {
          return true;
        }
      }
      if((seq[w].value[2] == seq[v].value[0]) && (seq[w].value[1] == seq[u].value[1]))
      {
        if(rule_4(seq, u, v, w, 2, 2))
        {
          return true;
        }
      }
      if((seq[w].value[2] == seq[v].value[0]) && (seq[w].value[1] == seq[u].value[2]))
      {
        if(rule_4(seq, u, v, w, 1, 2))
        {
          return true;
        }
      }
    }
  }
      
if((seq[v].value[2] == seq[u].value[0]) && (seq[v].value[1] != seq[u].value[1]) && (seq[v].value[1] != seq[u].value[2]))
  {
    if((seq[u].usd <= 1) || (seq[v].usd <= 1))
    {
      if((seq[w].value[1] == seq[v].value[0]) && (seq[w].value[2] == seq[u].value[1]))
      {
        if(rule_4(seq, u, v, w, 2, 1))
        {
          return true;
        }
      }
    
      if((seq[w].value[1] == seq[v].value[0]) && (seq[w].value[2] == seq[u].value[2]))
      {
        if(rule_4(seq, u, v, w, 1, 1))
        {
          return true;
        }
      }
      if((seq[w].value[2] == seq[v].value[0]) && (seq[w].value[1] == seq[u].value[1]))
      {
        if(rule_4(seq, u, v, w, 2, 1))
        {
          return true;
        }
      }
      if((seq[w].value[2] == seq[v].value[0]) && (seq[w].value[1] == seq[u].value[2]))
      {
        if(rule_4(seq, u, v, w, 1, 1))
        {
          return true;
        }
      }
    }
  }

  return false;   
}

bool rule_4(vector<list> &seq, int u, int v, int w, int idx_u, int idx_v)
{ 
  if((seq[u].usd <= 1) && (seq[u].flag < 0))
  { 
    if((seq[v].usd <= 1) && (seq[v].flag < 0))
    {
      seq[w].value[1] = seq[u].value[idx_u];
      seq[w].value[2] = seq[v].value[idx_v];
      seq.erase(seq.begin() + v);
      seq.erase(seq.begin() + u);
      return true;
    }
    else
    { 
      seq[w].value[1] = seq[u].value[idx_u]; 
      seq[w].value[2] = seq[v].value[idx_v];
      seq[v].value[1] = seq[w].value[0];

      if(idx_u == 1)
      {
		    seq[v].value[2] = seq[u].value[2];
      }
      else if(idx_u == 2)
      {
		    seq[v].value[2] = seq[u].value[1];
      }
     
      seq.insert(seq.begin() + v, seq[w]);//insert_list(seq, v, w);
      seq.erase(seq.begin() + w + 1);//delete_ele(seq, w + 1);
      seq.erase(seq.begin() +u);//delete_ele(seq, u);
      return true;
    }
  }
  else
    { 
      if((seq[v].usd <= 1) && (seq[v].flag < 0))
      {
        seq[w].value[1] = seq[u].value[idx_u];
        seq[w].value[2] = seq[v].value[idx_v];
        seq.erase(seq.begin() + v);//delete_ele(seq, v);
        return true;
      }
      return false;
    }
}
   
//t[u]==a+b, t[v]==c+d, t[w]==u+v   
bool reduce_rule_5(vector<list> &seq, int u, int v, int w)
{
  if((seq[w].value[1] == seq[u].value[0]) && (seq[w].value[2] == seq[v].value[0]))
  {
    if((seq[u].usd <= 1) || (seq[v].usd <= 1))   
    {
      if((seq[u].value [1] == seq[v].value[1]) && (seq[u].value[2] != seq[v].value[2]))
      {
        if(rule_5(seq, u, v, w, 2, 2))
        { 
            return true;
        }
      }
      if((seq[u].value[1] == seq[v].value[2]) && (seq[u].value[2] != seq[v].value[1]))
      {
        if(rule_5(seq, u, v, w, 2, 1))
        {
            return true;
        }  
      }

      if((seq[u].value[2] == seq[v].value[1]) && (seq[u].value[1] != seq[v].value[2]))
      {
        if(rule_5(seq, u, v, w, 1, 2))
        {
            return true;
        }
      }
      if((seq[u].value[2] == seq[v].value[2]) && (seq[u].value[1] != seq[v].value[1]))
      {
        if(rule_5(seq, u, v, w, 1, 1))
        {
            return true;
        }
      }
    }
  }

  if((seq[w].value [2] == seq[u].value[0]) && (seq[w].value[1] == seq[v].value[0]))
  {
    if((seq[u].usd <= 1) || (seq[v].usd <= 1))
    {
      if((seq[u].value[1] == seq[v].value[1]) && (seq[u].value[2] != seq[v].value[2]))
      {
        if(rule_5(seq, u, v, w, 2, 2))
        {
            return true;
        }
      }
      if((seq[u].value[1] == seq[v].value[2]) && (seq[u].value[2] != seq[v].value[1]))
      {
        if(rule_5(seq, u, v, w, 2, 1))
        {
            return true;
        }
      }
      if((seq[u].value[2] == seq[v].value[1]) && (seq[u].value[1] != seq[v].value[2]))
      {
        if(rule_5(seq, u, v, w, 1, 2))
        {
            return true;
        }
      }
      if((seq[u].value[2] == seq[v].value[2]) && (seq[u].value[1] != seq[v].value[1]))
      {
        if(rule_5(seq, u, v, w, 1, 1))
        {
            return true;
        }
      }
    }
  }

  return false;
}      
   
bool rule_5(vector<list> &seq, int u, int v, int w, int idx_u, int idx_v)
{ 
  if((seq[u].usd <= 1) && (seq[u].flag < 0))
  { 
    if((seq[v].usd <= 1) && (seq[v].flag < 0))
    {
      seq[w].value[1] = seq[u].value[idx_u];
      seq[w].value[2] = seq[v].value[idx_v];
      seq.erase(seq.begin() + v);
      seq.erase(seq.begin() + u);
      return true;
    }
    else
    { 
      seq[w].value[1] = seq[u].value[idx_u]; 
      seq[w].value[2] = seq[v].value[idx_v];
      seq.erase(seq.begin() + u);//delete_ele(seq, u);
      return true;
    }
  }
  else
    { 
      if((seq[v].usd <= 1) && (seq[v].flag < 0))
      {
        seq[w].value[1] = seq[u].value[idx_u];
        seq[w].value[2] = seq[v].value[idx_v];
        seq.erase(seq.begin() + v);//delete_ele(seq, v);
        return true;
      }
      return false;
    }
}


/*the 6 reduce_rule functions listed here have included the 11 rules shown in the paper*/
bool findreduce(vector<list> &seq)
{
  bool flag = false;
  for(int i = SIZE; i < seq.size(); i++)
  {
    for(int j = i + 1; j < seq.size(); j++)
    {
      if(redece_rule_1(seq, i, j)){return true;}
      for(int k = j + 1; k < seq.size(); k++)
      {
        if(reduce_rule_2(seq, i, j, k)){return true;}
        if(reduce_rule_3(seq, i, j, k)){return true;}
        if(reduce_rule_4(seq, i, j, k)){return true;}
        if(reduce_rule_5(seq, i, j, k)){return true;}
		if(reduce_rule_6(seq, i, j, k)){return true;}
      }
    }
  }
  
  return false;
}
