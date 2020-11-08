#include <vector>
#include "residue.h"
#include "LJPotential.h"
#include <Eigen/Dense>

using namespace std;

double calLJPotential(double d_star2, double eij, double lambd){
	double lj_potential = 0;
	if (d_star2 >= 0 && d_star2 <= 0.565){
		lj_potential = lambd*(49.69 - 40.06*sqrt(d_star2));
	}else if (d_star2 > 0.565 && d_star2 <= 0.797){
		lj_potential = lambd*eij*(pow(d_star2, -6) - 2 * pow(d_star2, -3));
	}else if (d_star2 > 0.797 && d_star2 < 6.25){
		lj_potential = eij*(pow(d_star2, -6) - 2 * pow(d_star2, -3));
	}else{
		lj_potential = 0;
	}

	return lj_potential;
}

double getAij2(Atom& a, Atom& b){
	int type1, type2;
	if(a.m_LJType > b.m_LJType){
		type1 = b.m_LJType;
		type2 = a.m_LJType;
	}else{
		type1 = a.m_LJType;
		type2 = b.m_LJType;	
	}
	double aij2 = -1;
	//Special cases for some atoms
	if(type1==9 || type1==10 || type1==11){
		if(type2==11 || type2==12){
			aij2=9.61;
		}else if(type2==13 || type2==14 || type2==15 || type2==16){
			aij2=8.41;
		}
	}else if(type1==12){
		if(type2==15 || type2==16){
			aij2=8.41;	
		}
	}else if(type1==13){
		if(type2==15 || type2==16){
			aij2=7.84;	
		}
	}else if(type1==15 || type1==16){
		if(type2==16){
			aij2=7.84;	
		}
	}else if(type1==17 && type2==17){
		aij2=4;
	}else if(type1==8 && type2==17){
		aij2=9;
	}

	if(aij2 == -1){
		aij2 = pow((a.m_radii + b.m_radii), 2);
	}

	return aij2;
}

double getMSLJPotential(vector<Residue>& residuesData, int res_id, int rl_id){

	double mslj_potentials = 0;
	//a atom side chain
	int rl_length = residuesData[res_id].m_rotamerlibs[rl_id].m_atoms.size();
	for (int i = 0; i < rl_length; i++){
		Vector3d* a_atom_sc = &residuesData[res_id].m_rotamerlibs[rl_id].m_atoms[i].m_position;
		//b atom main chain
		int cl_length = residuesData[res_id].m_contacts.size();
		for (int j = 0; j < cl_length; j++){
			int b_res_id = residuesData[res_id].m_contacts[j].m_j;
			//{ "N", "CA", "C", "O", "CB" }
			int end;
			if (residuesData[b_res_id].m_resname == 'G'){
				end = 4;
			}else{
				end = 5;
			}
			for (int k = 0; k < end; k++){
				Vector3d* b_atom_mc = &residuesData[b_res_id].m_atoms[LJ_MAIN[k]].m_position;
				//atom mightbe contact
				double dij2 = (*a_atom_sc - *b_atom_mc).squaredNorm();
				double aij2 = getAij2(residuesData[res_id].m_rotamerlibs[rl_id].m_atoms[i], residuesData[b_res_id].m_atoms[LJ_MAIN[k]]);
				double d_star2 = dij2 / aij2;
				double mslj_potential = 0;
				if (d_star2 < 6.25){
					double lambd = 1.6;
					if (residuesData[res_id].m_rotamerlibs[rl_id].m_atoms[i].m_LJType == 6 && residuesData[b_res_id].m_atoms[LJ_MAIN[k]].m_LJType == 6){
						lambd = 1;
					}
					double eij = sqrt(residuesData[res_id].m_rotamerlibs[rl_id].m_atoms[i].m_well*residuesData[b_res_id].m_atoms[LJ_MAIN[k]].m_well);
					mslj_potential = calLJPotential(d_star2, eij, lambd);
				}
				mslj_potentials = mslj_potentials + mslj_potential;
			}

		}
	
	}

	return mslj_potentials;

}

//double getSSLJPotential(vector<Residue>& residuesData, int res_id, int rl_id){
//
//	double sslj_potentials = 0;
//	//a atom side chain
//	int arl_length = residuesData[res_id].m_rotamerlibs[rl_id].m_atoms.size();
//	for (int i = 0; i < arl_length; i++){
//		Vector3d* a_atom_sc = &residuesData[res_id].m_rotamerlibs[rl_id].m_atoms[i].m_position;
//		//b atom side chain
//		int cl_length = residuesData[res_id].m_contacts.size();
//		for (int j = 0; j < cl_length; j++){
//			int b_res_id = residuesData[res_id].m_contacts[j].m_j;
//			if (residuesData[b_res_id].m_resname == 'G' || residuesData[b_res_id].m_resname == 'A'){
//				continue;
//			}
//			int b_rl_id = residuesData[b_res_id].m_rl_id;
//			int brl_length = residuesData[b_res_id].m_rotamerlibs[b_rl_id].m_atoms.size();
//			for (int k = 0; k < brl_length; k++){
//				Vector3d* b_atom_mc = &residuesData[b_res_id].m_rotamerlibs[b_rl_id].m_atoms[k].m_position;
//				//atom mightbe contact
//				double dij2 = (*a_atom_sc - *b_atom_mc).squaredNorm();
//				double aij2 = getAij2(residuesData[res_id].m_rotamerlibs[rl_id].m_atoms[i], residuesData[b_res_id].m_rotamerlibs[b_rl_id].m_atoms[k]);
//				double d_star2 = dij2 / aij2;
//				double sslj_potential = 0;
//				if (d_star2 < 6.25){
//					double lambd = 1.6;
//					if (residuesData[res_id].m_rotamerlibs[rl_id].m_atoms[i].m_LJType == 6 && residuesData[b_res_id].m_rotamerlibs[b_rl_id].m_atoms[k].m_LJType == 6){
//						lambd = 1;
//					}
//					double eij = sqrt(residuesData[res_id].m_rotamerlibs[rl_id].m_atoms[i].m_well*residuesData[b_res_id].m_rotamerlibs[b_rl_id].m_atoms[k].m_well);
//					sslj_potential = calLJPotential(d_star2, eij, lambd);
//				}
//				sslj_potentials = sslj_potentials + sslj_potential;
//			}
//
//		}
//
//	}
//
//	return sslj_potentials;
//
//}

double getSSLJPotential(vector<Residue>& residuesData, int a_res_id, int a_rl_id){

	double sslj_potentials = 0;
	int cm_length = residuesData[a_res_id].m_contacts.size();
	for (int i = 0; i < cm_length; i++){
		double sslj_rr_potentials = 0;
		int b_res_id = residuesData[a_res_id].m_contacts[i].m_j;
		int b_rl_id = residuesData[b_res_id].m_rl_id;
		if (residuesData[b_res_id].m_resname == 'G' || residuesData[b_res_id].m_resname == 'A'){
			continue;
		}
		int a_key = b_res_id *1000 + b_rl_id;
		//load cash
		if (residuesData[a_res_id].m_rotamerlibs[a_rl_id].m_rr_sslj.count(a_key)>0){
			sslj_potentials = sslj_potentials + residuesData[a_res_id].m_rotamerlibs[a_rl_id].m_rr_sslj[a_key];
			continue;
		}
		//a atom side chain
		int a_atoms_length = residuesData[a_res_id].m_rotamerlibs[a_rl_id].m_atoms.size();
		for (int j = 0; j < a_atoms_length; j++){
			Vector3d* a_atom_sc = &residuesData[a_res_id].m_rotamerlibs[a_rl_id].m_atoms[j].m_position;
			//b atom side chain
			int b_atoms_length = residuesData[b_res_id].m_rotamerlibs[b_rl_id].m_atoms.size();
			for (int k = 0; k < b_atoms_length; k++){
				Vector3d* b_atom_sc = &residuesData[b_res_id].m_rotamerlibs[b_rl_id].m_atoms[k].m_position;
				//atom mightbe contact
				double dij2 = (*a_atom_sc - *b_atom_sc).squaredNorm();
				double aij2 = getAij2(residuesData[a_res_id].m_rotamerlibs[a_rl_id].m_atoms[j], residuesData[b_res_id].m_rotamerlibs[b_rl_id].m_atoms[k]);
				double d_star2 = dij2 / aij2;
				double sslj_potential = 0;
				if (d_star2 < 6.25){
					double lambd = 1.6;
					if (residuesData[a_res_id].m_rotamerlibs[a_rl_id].m_atoms[j].m_LJType == 6 && residuesData[b_res_id].m_rotamerlibs[b_rl_id].m_atoms[k].m_LJType == 6){
						lambd = 1;
					}
					double eij = sqrt(residuesData[a_res_id].m_rotamerlibs[a_rl_id].m_atoms[j].m_well*residuesData[b_res_id].m_rotamerlibs[b_rl_id].m_atoms[k].m_well);
					sslj_potential = calLJPotential(d_star2, eij, lambd);
				}
				//sslj_potentials = sslj_potentials + sslj_potential;
				sslj_rr_potentials = sslj_rr_potentials + sslj_potential;
			}
		}
		sslj_potentials = sslj_potentials + sslj_rr_potentials;
		//save to cash
		residuesData[a_res_id].m_rotamerlibs[a_rl_id].m_rr_sslj[a_key] = sslj_rr_potentials;
		int b_key = a_res_id*1000 + a_rl_id;
		residuesData[b_res_id].m_rotamerlibs[b_rl_id].m_rr_sslj[b_key] = sslj_rr_potentials;
	}
	return sslj_potentials;
}

