#ifndef OTHERRL_H
#define OTHERRL_H

#include <string> 
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <vector>
#include <sstream>
#include <fstream> 
#include <iostream> 
#include <map>

#include "residue.h"
#include "config.h"
#include "pdb.h"
#include "angle.h"

using namespace std;

map<string, string> getConstraintList(string& path);

void addRotamersInfoFromFASPR(const string& inputname, const string& tmp_dir, vector<Residue>& residuesData);

void addRotamersInfoFromOthers(const string& inputname, const map<string, string>& other_list, vector<Residue>& residuesData);

void addRotamersInfoFromOSCAR(string inputname, vector<Residue>& residuesData);

void addRotamersInfoFromRotaNN(const string& inputname, const map<string, string>& rotann_list, vector<Residue>& residuesData);

void calOSCAR(const string& list_path);

void addRotamersInfoRefineX4(vector<Residue>& residuesData);

#endif