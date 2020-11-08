#include <iostream> 
#include<string> 
#include<vector> 
#include <map>
#include "atom.h"
#include "residue.h"

using namespace std;
using namespace Eigen;

DASFReference::DASFReference(){
}



Vector3d DASFReference::transCoordinate(Vector3d& position){
	Vector3d position_new = position - ref;
	Vector4d position_new4d(position_new[0], position_new[1], position_new[2], 1);
	Vector4d r = position_new4d.transpose()*rotation_matrix;
	return Vector3d(r[0],r[1],r[2]);
}

Contact::Contact(int i, int j, double bb_distance){
	m_i = i;
	m_j = j;
	m_bb_distance = bb_distance;
}

RotamerLib::RotamerLib(int id, double prob, double x1, double x2, double x3, double x4){
	m_id = id;
	m_prob = prob;
	m_dihedral.push_back(x1);
	m_dihedral.push_back(x2);
	m_dihedral.push_back(x3);
	m_dihedral.push_back(x4);
}


DASFFeature::DASFFeature(vector<double> mean, vector<double> sd, int nums, int window_len){
	m_mean = Vector3d(mean[0], mean[1], mean[2]);
	m_sd = Vector3d(sd[0], sd[1], sd[2]);
	m_nums = nums;
	m_window_len = window_len;
}

Residue::Residue(){
	m_rl_id = 0;
}

Residue::Residue(int param_resid, char param_resname){
	m_rl_id = 0;
	m_resid = param_resid;
	m_resname = param_resname;
}

void Residue::getDASFReference(){
	Vector3d ca_ref = m_atoms["CA"].m_position;
	Vector3d c_ref = m_atoms["C"].m_position;
	Vector3d o_ref = m_atoms["O"].m_position;

	m_DASFReference.ref = ca_ref;
	Vector3d c_ref_new = c_ref - m_DASFReference.ref;
	Vector3d o_ref_new = o_ref - m_DASFReference.ref;

	//c - ca
	Vector3d x_axis = c_ref_new / c_ref_new.norm();
	Vector3d c_o = o_ref_new - c_ref_new;

	//o - c perpendicular to x_axis
	Vector3d y_axis = c_o - (x_axis.dot(c_o) / x_axis.dot(x_axis) * x_axis);
	y_axis = y_axis / y_axis.norm();

	Vector3d z_axis = x_axis.cross(y_axis);

	m_DASFReference.rotation_matrix << x_axis[0], y_axis[0], z_axis[0], 0, 
									   x_axis[1], y_axis[1], z_axis[1], 0,
									   x_axis[2], y_axis[2], z_axis[2], 0,
									   0, 0, 0, 1;

}

Atom Residue::getAtom(const string& atomname){
	if(m_atoms.count(atomname)>0){
		return m_atoms[atomname];
	}else{
		throw atomname + "lost at: " + to_string(m_resid) + " " + m_resname;
	}
}

void Residue::setDihedrals(vector<double> param_dihedrals){
	m_dihedrals = param_dihedrals;
}

string getTriResname(char resname){
	string resnameTri;
	switch(resname){
		case('G'):
			resnameTri = "GLY";
			break;
		case('A'):
			resnameTri = "ALA";
			break;
		case('S'):
			resnameTri = "SER";
			break;
		case('C'):
			resnameTri = "CYS";
			break;
		case('V'):
			resnameTri = "VAL";
			break;
		case('I'):
			resnameTri = "ILE";
			break;
		case('L'):
			resnameTri = "LEU";
			break;
		case('T'):
			resnameTri = "THR";
			break;
		case('R'):
			resnameTri = "ARG";
			break;
		case('K'):
			resnameTri = "LYS";
			break;
		case('D'):
			resnameTri = "ASP";
			break;
		case('E'):
			resnameTri = "GLU";
			break;
		case('N'):
			resnameTri = "ASN";
			break;
		case('Q'):
			resnameTri = "GLN";
			break;
		case('M'):
			resnameTri = "MET";
			break;
		case('H'):
			resnameTri = "HIS";
			break;
		case('P'):
			resnameTri = "PRO";
			break;
		case('F'):
			resnameTri = "PHE";
			break;		
		case('Y'):
			resnameTri = "TYR";
			break;
		case('W'):
			resnameTri = "TRP";
			break;
		default:
			cout << resname << endl;
			throw "residue.getTriResname() wrong";
	}
	return resnameTri;
}

char getResname(const string& resnameTri){
	char resname;
	if(resnameTri.find("GLY") != string::npos){
		resname = 'G';
	}else if(resnameTri.find("ALA") != string::npos){
		resname = 'A';
	}else if(resnameTri.find("SER") != string::npos){
		resname = 'S';
	}else if(resnameTri.find("CYS") != string::npos){
		resname = 'C';
	}else if(resnameTri.find("VAL") != string::npos){
		resname = 'V';
	}else if(resnameTri.find("ILE") != string::npos){
		resname = 'I';
	}else if(resnameTri.find("LEU") != string::npos){
		resname = 'L';
	}else if(resnameTri.find("THR") != string::npos){
		resname = 'T';
	}else if(resnameTri.find("ARG") != string::npos){
		resname = 'R';
	}else if(resnameTri.find("LYS") != string::npos){
		resname = 'K';
	}else if(resnameTri.find("ASP") != string::npos){
		resname = 'D';
	}else if(resnameTri.find("GLU") != string::npos){
		resname = 'E';
	}else if(resnameTri.find("ASN") != string::npos){
		resname = 'N';
	}else if(resnameTri.find("GLN") != string::npos){
		resname = 'Q';
	}else if(resnameTri.find("MET") != string::npos){
		resname = 'M';
	}else if(resnameTri.find("HIS") != string::npos){
		resname = 'H';
	}else if(resnameTri.find("PRO") != string::npos){
		resname = 'P';
	}else if(resnameTri.find("PHE") != string::npos){
		resname = 'F';
	}else if(resnameTri.find("TYR") != string::npos){
		resname = 'Y';
	}else if(resnameTri.find("TRP") != string::npos){
		resname = 'W';
	}else{

		cout << resnameTri << endl;
		throw "residue.getResname() wrong";
	}
	return resname;

	
}

void addAtom(vector<Residue>& residuesData, const Atom& atom){

	int current_length = residuesData.size();
	int resid = atom.m_resid;
	char resname = atom.m_resname;

	//first residue
	if(current_length == 0){
		Residue r(resid, resname);
		r.m_atoms.insert(pair<string, Atom>(atom.m_name,atom));
		residuesData.push_back(r);
	}else{
		if(resid == residuesData[current_length-1].m_resid){
			residuesData[current_length-1].m_atoms.insert(pair<string, Atom>(atom.m_name,atom));
		}else{
			Residue r(resid, resname);
			r.m_atoms.insert(pair<string, Atom>(atom.m_name,atom));
			residuesData.push_back(r);		
		}
	} 
}

vector<Residue> atomToResidue(const vector<Atom>& atomsData, bool main_chain){
	vector<Residue> residuesData;
	int length = atomsData.size();
	for(int i=0; i<length; i++){
		char restype = atomsData[i].m_restype;
		if(main_chain && !atomsData[i].m_isMainChain){
			continue;
		}
		if(atomsData[i].m_LJType == -1 || restype != 'A'){
			//cout << "No use atom at: " + to_string(atomsData[i].m_resid) + " " + atomsData[i].m_resname + " " + atomsData[i].m_name  << endl;
			continue;
		}
		addAtom(residuesData, atomsData[i]);
	}
	return residuesData;

}
