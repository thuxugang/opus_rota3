#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector> 
#include <iostream>
#include <sstream>
#include "residue.h"
#include "redis.h"
#include "geometry.h"
#include "hiredis/hiredis.h"

using namespace std;

string DB_ROTAMER;
vector<string> DB_DASFS;
vector<int> DASFS_KEYLEN;

vector<string> splits(const string &s, char delim){
	stringstream ss(s);
	string item;
	vector<string> elems;
	while (getline(ss, item, delim)){
		elems.push_back(item);
	}
	return elems;
}

vector<double> splitd(const string &s, char delim){
	stringstream ss(s);
	string item;
	vector<double> elems;
	while (getline(ss, item, delim)){
		elems.push_back(atof(item.c_str()));
	}
	return elems;
}

Redis::Redis(){

	c = redisConnect(CONFIG["ip"].c_str(), stoi(CONFIG["port"]));
	if (c->err) {
		cout << "Failed to connect redis!" << endl;
		throw "Failed to connect redis!";
	}
	DB_ROTAMER = CONFIG["db_rotamer"];
	if (CONFIG["db_dasf_5"] != ""){
		DB_DASFS.push_back(CONFIG["db_dasf_5"]);
		DASFS_KEYLEN.push_back(5);
	}
	if (CONFIG["db_dasf_7"] != ""){
		DB_DASFS.push_back(CONFIG["db_dasf_7"]);
		DASFS_KEYLEN.push_back(7);
	}
	if (CONFIG["db_dasf_9"] != ""){
		DB_DASFS.push_back(CONFIG["db_dasf_9"]);
		DASFS_KEYLEN.push_back(9);
	}
	if (CONFIG["db_dasf_11"] != ""){
		DB_DASFS.push_back(CONFIG["db_dasf_11"]);
		DASFS_KEYLEN.push_back(11);
	}
}

void Redis::destroy(){

	redisFree(c);
}

void Redis::addAllRotamers(Residue& residue, int top_n){

	try{
		redisReply* reply = (redisReply *)redisCommand(c, "GET %s", residue.m_phipsi.c_str());
		if (reply->type != REDIS_REPLY_NIL){
			string result = reply->str;
			vector<string> rls = splits(result, ';');
			int length_rl = rls.size();
			int max_len = top_n < length_rl ? top_n : length_rl;
			for (int i = 0; i < max_len; i++){
				vector<double> rotamers = splitd(rls[i], '_');
				assert(rotamers.size() == 5);
				RotamerLib rl(i, rotamers[0], rotamers[1], rotamers[2], rotamers[3], rotamers[4]);
				residue.m_rotamerlibs.push_back(rl);
			}
		
		}
		freeReplyObject(reply);

	}catch(...){
		cout << "rotamer query failed!" << endl;
		throw "rotamer query failed!";
	}
}

int Redis::addAllDASFs(string& key, vector<Residue>& residuesData, vector<Geo*>& geosData, int center, int half, int window_len){

	vector<DASFFeature> fs;
	try{
		redisReply* reply = (redisReply *)redisCommand(c, "GET %s", key.c_str());
		if (reply->type != REDIS_REPLY_NIL){
			string result = reply->str;
			vector<string> features = splits(result, ';');
			int length_features = features.size();
			for (int i = 0; i < length_features; i++){
				vector<string> feature = splits(features[i], '#');
				int length_feature = feature.size();
				fs.push_back(DASFFeature(splitd(feature[0], '_'), splitd(feature[1], '_'), atoi(feature[2].c_str()), window_len));
			}
		}
		freeReplyObject(reply);
	}
	catch (...){
		cout << "dasf query failed!" << endl;
		throw "dasf query failed!";
	}

	if(fs.size() != 0){
		int start;
		if(half%2 ==0){
			start = 0;
		}else{
			start = 1;
		}
		int f_id = 0;
		for(int i = center-half+start;i <= center+half-start; i=i+2){
			//init
			if (residuesData[i].m_standerDASFs.size() == 0){
				for (int j = 0; j < geosData[i]->m_num_dasfs; j++){
					vector<DASFFeature> fs;
					residuesData[i].m_standerDASFs.push_back(fs);
				}
			}
			for(int j=0; j < geosData[i]->m_num_dasfs;j++){
				residuesData[i].m_standerDASFs[j].push_back(fs[f_id]);
				f_id = f_id + 1;		
			}				
		}
		assert(f_id == fs.size());	
	}
	if (fs.size() != 0){
		return 1;
	}else{
		return 0;
	}
}

void addRotamersInfo(vector<Residue>& residuesData, Redis& r, int top_n){
	r.m_rltable = DB_ROTAMER;
	redisReply *pRedisReply = (redisReply *)redisCommand(r.c, "SELECT %s", r.m_rltable.c_str());
	int length = residuesData.size();
	for(int i=0; i<length; i++){
		r.addAllRotamers(residuesData[i], top_n);
	}
}


float addDASFInfo(vector<Residue>& residuesData, vector<Geo*>& geosData, Redis& r){

	int length = DB_DASFS.size();
	float totals = 0;
	float contains = 0;
	for (int i = 0; i<length; i++){
		r.m_rltable = DB_DASFS[i];
		(redisReply *)redisCommand(r.c, "SELECT %s", r.m_rltable.c_str());
		int window_len = DASFS_KEYLEN[i];
		int num = 0;
		int half_1 = int((window_len-1)/2);
		int half_2 = int((window_len+1)/2);
		int length_seq = residuesData.size();
		for (int j = half_1; j< length_seq-half_1; j++){
			string key = "";
			int m = half_1;
			while (m>0){
				key = key + residuesData[j-m].m_resname;
				m = m - 1;
			}
			key = key + residuesData[j].m_resname;
			int n = 1;
			while (n<half_2){
				key = key + residuesData[j+n].m_resname;
				n = n + 1;
			}
			totals += 1;
			contains += r.addAllDASFs(key, residuesData, geosData, j, half_1, window_len);
		}	
	}
	return contains / totals;
}
