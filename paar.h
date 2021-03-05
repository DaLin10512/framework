#ifndef _PAAR_H_
#define _PAAR_H_

#include"const.h"
#include"regidster.h"
#include"bfi.h"
vector<pair<int, int>> find_max_idx(uint64_t *input_matrix, int number_of_columns);
int Paar_algorithm1(vector<list> &seq, uint64_t NumInputs, uint64_t NumTargets, uint64_t *Target);
int RPaar_algorithm1(vector<list> &seq, uint64_t NumInputs, uint64_t NumTargets, uint64_t *Target);

#endif
