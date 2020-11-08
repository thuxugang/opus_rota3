#ifndef REDIS_H
#define REDIS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector> 
#include "residue.h"
#include "geometry.h"
#include "config.h"
#include "hiredis/hiredis.h"

using namespace std;

class Redis{
public:

	string m_rltable;
	redisContext *c;

	Redis();
	void addAllRotamers(Residue& residue, int top_n);
	int addAllDASFs(string& key, vector<Residue>& residuesData, vector<Geo*>& geosData, int center, int half, int window_len);
	void destroy();
};


void addRotamersInfo(vector<Residue>& residuesData, Redis& r, int top_n);
float addDASFInfo(vector<Residue>& residuesData, vector<Geo*>& geosData, Redis& r);



#endif