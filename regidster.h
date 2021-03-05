
#ifndef _REGIDSTER_H_
#define _REGIDSTER_H_
#include"const.h"

typedef struct 
{
  vector<int> value;
  int flag;
  int usd;
}list;

typedef struct
{
    vector<list> seq;
    int gap;
    int start;
    //int len;
}thread_data;

bool findreduce(vector<list> &seq);
void rule_1(vector<list> &seq, int u, int v, int ind_u);
void rule_2(vector<list> &seq, int u, int v, int w, int idx_u);
void rule_3(vector<list> &seq, int u, int v, int w, int idx_v);
bool rule_4(vector<list> &seq, int u, int v, int w, int idx_u, int idx_v);
bool rule_5(vector<list> &seq, int u, int v, int w, int idx_u, int idx_v);
void rule_6(vector<list> &seq, int u, int v, int w, int idx_u);
bool find_multi(vector<list> &seq);
void update_seq(vector<list> &seq);
bool delete_twice_ele(vector<int> &aray);
vector<uint64_t> find_targets(vector<list> seq);
void sort_seq(vector<list> &seq);
bool pre_reduce_seq(vector<list> &seq);
void pre_bp(vector<list> &seq, int sp, int gap);
void count_used_time(vector<list> &seq);


#endif




