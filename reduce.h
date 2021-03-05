#ifndef _READ_H_
#define _READ_H_
#include"const.h"
#include"regidster.h"


bool check(vector<list> seq);
void write_file(vector<list> seq);
int get_reduced(vector<list> &seq);
vector<list> Algorithm_B(vector<list> e_seq, vector<int> input, uint64_t *Target, uint64_t *Dist, uint64_t *InitDist);
void get_equivalent_circuit(vector<list> &seq, int sp, int length);
void find_inputs(vector<list> tmp_seq, vector<int> &input);
void find_outputs(vector<list> tmp_seq, vector<int> input, uint64_t *Target, uint64_t *Dist, uint64_t *InitDist);

#endif



