#include <iostream> 
#include <random>
#include <string> 
#include <vector> 
#include "residue.h"
#include "geometry.h"
#include "MonteCarlo.h"

using namespace std;


void MC_init(vector<Residue>& residuesData){
	int length = residuesData.size();
	for(int i=0;i<length; i++){
		if(residuesData[i].m_resname != 'G' && residuesData[i].m_resname != 'A'){
			int rl_length = residuesData[i].m_rotamerlibs.size();
			int min_id = 0;
			double min_score = residuesData[i].m_rotamerlibs[0].m_totalScores;
			for(int j=0;j<rl_length;j++){
				if(residuesData[i].m_rotamerlibs[j].m_totalScores < min_score){
					min_id = j;
					min_score = residuesData[i].m_rotamerlibs[j].m_totalScores;
				}
			}
			residuesData[i].m_rl_id = min_id;
			//cout << residuesData[i].m_resid << " " << residuesData[i].m_resname << " " <<  min_id << " " << 
			//	residuesData[i].m_rotamerlibs[min_id].m_DASFScores << " " << residuesData[i].m_rotamerlibs[min_id].m_RLScores
			//	<< " " << residuesData[i].m_rotamerlibs[min_id].m_MS_LJScores << " " << residuesData[i].m_rotamerlibs[min_id].m_totalScores << endl;
		}
		
	}

}

void MC_simulation(vector<Residue>& residuesData){
	int res_length = residuesData.size();
	vector<int> idx;
	for (int i = 0; i < res_length; i++){
		idx.push_back(i);
	}
	int optimization_times = stoi(CONFIG["optimization_times"]);
	for (int i = 0; i < optimization_times; i++){
		shuffle(idx.begin(), idx.end(), E);
		for (int j = 0; j < res_length; j++){
			int id = idx[j];
			//cout << id << residuesData[id].m_resname << endl;
			if (residuesData[id].m_resname != 'G' && residuesData[id].m_resname != 'A'){
				int rl_length = residuesData[id].m_rotamerlibs.size();
				int rl_old = residuesData[id].m_rl_id;
				uniform_int_distribution<unsigned> u(0, rl_length - 1);
				int rl_new = u(E);
				//exclude same
				if (rl_old != rl_new){
					double sslj_old = stod(CONFIG["w_sslj"])*getSSLJPotential(residuesData, id, rl_old);
					double score_old = residuesData[id].m_rotamerlibs[rl_old].m_totalScores + sslj_old;
					double sslj_new = stod(CONFIG["w_sslj"])*getSSLJPotential(residuesData, id, rl_new);
					double score_new = residuesData[id].m_rotamerlibs[rl_new].m_totalScores + sslj_new;
					if (score_new < score_old){
						residuesData[id].m_rl_id = rl_new;
					}
				}
			}
		}
	}
}

void getBestStructure(vector<Residue>& residuesData){
	//get lowest m_totalScores
	MC_init(residuesData);

	//simulation
	if(stoi(CONFIG["optimized"])){
		MC_simulation(residuesData);
	}
	
}