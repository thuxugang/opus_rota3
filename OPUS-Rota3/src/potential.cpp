#include <iostream> 
#include <vector> 
#include "residue.h"
#include "geometry.h"
#include "potential.h"
#include "rotamerLibScore.h"
#include "LJPotential.h"
#include "DASFScore.h"
#include "config.h"

using namespace std;

vector<Contact> getContactList(const vector<Geo*>& geosData, vector<Residue>& residuesData){
	vector<Contact> contact_list;
	int seq_length = geosData.size();
	for(int i=0; i<seq_length; i++){
		Vector3d* a_cb;
		if(geosData[i]->m_resname == 'G'){
			a_cb = &geosData[i]->CA.m_position;
		}else{
			a_cb = &geosData[i]->CB.m_position;
		}
		//include self
		for(int j=i+1; j<seq_length; j++){
			Vector3d* b_cb;
			if(geosData[j]->m_resname == 'G'){
				b_cb = &geosData[j]->CA.m_position;
			}else{
				b_cb = &geosData[j]->CB.m_position;
			}		
			double bb_distance = (*a_cb-*b_cb).norm();
			double max_distance = bb_distance - geosData[i]->m_maxdis - geosData[j]->m_maxdis;
			//2.5*(2+2)
			if(max_distance < 10){
				contact_list.push_back(Contact(i,j,bb_distance));
				residuesData[i].m_contacts.push_back(Contact(i,j,bb_distance));
				residuesData[j].m_contacts.push_back(Contact(j,i,bb_distance));
			}
		}
	
	}
	return contact_list;

}


void initAllPotentials(vector<Geo*>& geosData, vector<Residue>& residuesData, double coverage){
	assert(geosData.size() == residuesData.size());
	int seq_length = geosData.size();
	for(int i=0; i<seq_length; i++){
		//cout << residuesData[i].m_resid << residuesData[i].m_resname << endl;
		int rl_length = residuesData[i].m_rotamerlibs.size();
		residuesData[i].getDASFReference();
		for(int j=0; j<rl_length; j++){
			RotamerLib* rl = &residuesData[i].m_rotamerlibs[j];
			geosData[i]->setRotamers(rl->m_dihedral);
			geosData[i]->addSideChain();
			geosData[i]->storeSideChain(*rl);

			//rotamer lib score
			residuesData[i].m_rotamerlibs[j].m_RLScores = stod(CONFIG["w_rp"])*getRLScore(residuesData[i], j) + \
				stod(CONFIG["w_x1cons"])*abs(residuesData[i].m_rotamerlibs[j].m_dihedral[0] - residuesData[i].m_rotamerlibs[0].m_dihedral[0]) + \
				stod(CONFIG["w_x2cons"])*abs(residuesData[i].m_rotamerlibs[j].m_dihedral[1] - residuesData[i].m_rotamerlibs[0].m_dihedral[1]);
			//main chain - side chain LJ potential
			residuesData[i].m_rotamerlibs[j].m_MS_LJScores = stod(CONFIG["w_mslj"])*getMSLJPotential(residuesData, i, j);
			//DASF score
			residuesData[i].m_rotamerlibs[j].m_DASFScores = coverage*stod(CONFIG["w_dasf"])*getDASFPotential(residuesData, i, j);
			//total
			residuesData[i].m_rotamerlibs[j].m_totalScores = residuesData[i].m_rotamerlibs[j].m_RLScores+residuesData[i].m_rotamerlibs[j].m_MS_LJScores+residuesData[i].m_rotamerlibs[j].m_DASFScores;			
		}
	}

}

