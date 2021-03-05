#ifndef _RANPER_H_
#define _RANPER_H_

#include"const.h"
#include"regidster.h"
vector<uint64_t> generate_p(int number);
vector<uint64_t> p_m_q(uint64_t *Target, uint64_t no_row, uint64_t no_col, uint64_t *Dist, vector<uint64_t> P, vector<uint64_t> Q);
vector<uint64_t> p_m_q_paar(uint64_t *Target, uint64_t no_row, uint64_t no_col, vector<uint64_t> P, vector<uint64_t> Q);
void recover_implementation(vector<list> &seq, vector<uint64_t> P, vector<uint64_t> Q);
#endif
