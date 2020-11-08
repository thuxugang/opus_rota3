#include <vector> 
#include "residue.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

vector<Vector3d> getFeatures(Residue& residue, int rl_id){
	vector<Vector3d> fs;
	switch (residue.m_resname){
		case 'V':
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[0].m_position);	//"CG1"
			break;
		case 'I':
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[0].m_position);	//"CG1"
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[2].m_position);	//"CD"
			break;
		case 'L':
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[0].m_position);	//"CG"
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[1].m_position);	//"CD1"
			break;
		case 'S':
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[0].m_position);	//"OG"
			break;
		case 'T':
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[0].m_position);	//"OG1"
			break;                            
		case 'D':
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[0].m_position);	//"CG"
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[1].m_position);	//"OD1"
			break;
		case 'N':
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[0].m_position);	//"CG"
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[1].m_position);	//"OD1"
			break;
		case 'E':
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[0].m_position);	//"CG"
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[1].m_position);	//"CD"
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[2].m_position);	//"OE1"
			break;
		case 'Q':
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[0].m_position);	//"CG"
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[1].m_position);	//"CD"
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[2].m_position);	//"OE1"
			break;
		case 'K':
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[0].m_position);	//"CG"
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[1].m_position);	//"CD"
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[2].m_position);	//"CE"
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[3].m_position);	//"NZ"
			break;
		case 'R':
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[0].m_position);	//"CG"
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[1].m_position);	//"CD"
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[2].m_position);	//"NE"
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[3].m_position);	//"CZ"
			break;
		case 'C':
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[0].m_position);	//"SG"
			break;
		case 'M':
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[0].m_position);	//"CG"
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[1].m_position);	//"SD"
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[2].m_position);	//"CE"
			break;
		case 'F':
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[0].m_position);	//"CG"
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[1].m_position);	//"CD1"
			break;
		case 'Y':
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[0].m_position);	//"CG"
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[1].m_position);	//"CD1"
			break;
		case 'W':
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[0].m_position);	//"CG"
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[1].m_position);	//"CD1"
			break;
		case 'H':
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[0].m_position);	//"CG"
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[1].m_position);	//"ND1"
			break;
		case 'P':
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[0].m_position);	//"CG"
			fs.push_back(residue.m_rotamerlibs[rl_id].m_atoms[1].m_position);	//"CD"
			break;
	}
	return fs;
}

double getDASFPotential(vector<Residue>& residuesData, int res_id, int rl_id){
	double results = 0;
	if (residuesData[res_id].m_resname != 'G' && residuesData[res_id].m_resname != 'A'){
		if (residuesData[res_id].m_standerDASFs.size() == 0){
			return results;
		}
		vector<Vector3d> fs = getFeatures(residuesData[res_id], rl_id);
		int fs_length = fs.size();
		for (int i = 0; i < fs_length; i++){
			Vector3d f = residuesData[res_id].m_DASFReference.transCoordinate(fs[i]);
			vector<DASFFeature> fs_standers = residuesData[res_id].m_standerDASFs[i];
			int fs_stander_length = fs_standers.size();
			for (int j = 0; j < fs_stander_length; j++){
				Vector3d v = f - fs_standers[j].m_mean;
				double result = fs_standers[j].m_window_len*(abs(v[0] / fs_standers[j].m_sd[0]) + abs(v[1] / fs_standers[j].m_sd[1]) + abs(v[2] / fs_standers[j].m_sd[2]))/10;
				results = results + result;
			}
		}
	}
	return results;
}