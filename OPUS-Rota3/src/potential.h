#ifndef POTENTIAL_H
#define POTENTIAL_H

#include<vector> 
#include "residue.h"
#include "geometry.h"
#include "config.h"

using namespace std;

vector<Contact> getContactList(const vector<Geo*>& geosData, vector<Residue>& residuesData);

void initAllPotentials(vector<Geo*>& geosData, vector<Residue>& residuesData, double coverage);

#endif