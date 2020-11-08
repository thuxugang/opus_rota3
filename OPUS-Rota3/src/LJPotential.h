#ifndef LJPOTENTIAL_H
#define LJPOTENTIAL_H

#include <vector> 
#include "residue.h"

double getMSLJPotential(vector<Residue>& residuesData, int res_id, int rl_id);
double getSSLJPotential(vector<Residue>& residuesData, int res_id, int rl_id);

#endif