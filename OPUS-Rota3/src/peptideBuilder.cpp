#include <iostream> 
#include <Eigen/Dense>
#include <string>
#include <vector> 
#include "angle.h"
#include "residue.h"
#include "atom.h"
#include "geometry.h"
#include "peptideBuilder.h"

using namespace std;
using namespace Eigen;


vector<Residue> initResidueChain(const vector<char>& resnames){
	vector<Residue> residuesData;
	int length = resnames.size();
	for(int i=0; i<length; i++){
		Residue aa = Residue(i+1,resnames[i]);
		residuesData.push_back(aa);
	}
	return residuesData;
}

Geo* getGeo(int resid, char resname, int* totalAtoms){
	Geo* geo;
	switch(resname){
	case 'G':
		geo = new GlyGeo(resid, *totalAtoms);
		*totalAtoms = *totalAtoms + 4;
		break;
	case 'A':
		geo = new AlaGeo(resid, *totalAtoms);
		*totalAtoms = *totalAtoms + 5;
		break;
	case 'S':
		geo = new SerGeo(resid, *totalAtoms);
		*totalAtoms = *totalAtoms + 6;
		break;
	case 'C':
		geo = new CysGeo(resid, *totalAtoms);
		*totalAtoms = *totalAtoms + 6;
		break;
	case 'V':
		geo = new ValGeo(resid, *totalAtoms);
		*totalAtoms = *totalAtoms + 7;
		break;
	case 'I':
		geo = new IleGeo(resid, *totalAtoms);
		*totalAtoms = *totalAtoms + 8;
		break;
	case 'L':
		geo = new LeuGeo(resid, *totalAtoms);
		*totalAtoms = *totalAtoms + 8;
		break;
	case 'T':
		geo = new ThrGeo(resid, *totalAtoms);
		*totalAtoms = *totalAtoms + 7;
		break;
	case 'R':
		geo = new ArgGeo(resid, *totalAtoms);
		*totalAtoms = *totalAtoms + 11;
		break;
	case 'K':
		geo = new LysGeo(resid, *totalAtoms);
		*totalAtoms = *totalAtoms + 9;
		break;
	case 'D':
		geo = new AspGeo(resid, *totalAtoms);
		*totalAtoms = *totalAtoms + 8;
		break;
	case 'N':
		geo = new AsnGeo(resid, *totalAtoms);
		*totalAtoms = *totalAtoms + 8;
		break;
	case 'E':
		geo = new GluGeo(resid, *totalAtoms);
		*totalAtoms = *totalAtoms + 9;
		break;
	case 'Q':
		geo = new GlnGeo(resid, *totalAtoms);
		*totalAtoms = *totalAtoms + 9;
		break;
	case 'M':
		geo = new MetGeo(resid, *totalAtoms);
		*totalAtoms = *totalAtoms + 8;
		break;
	case 'H':
		geo = new HisGeo(resid, *totalAtoms);
		*totalAtoms = *totalAtoms + 10;
		break;
	case 'P':
		geo = new ProGeo(resid, *totalAtoms);
		*totalAtoms = *totalAtoms + 7;
		break;
	case 'F':
		geo = new PheGeo(resid, *totalAtoms);
		*totalAtoms = *totalAtoms + 11;
		break;
	case 'Y':
		geo = new TyrGeo(resid, *totalAtoms);
		*totalAtoms = *totalAtoms + 12;
		break;
	case 'W':
		geo = new TrpGeo(resid, *totalAtoms);
		*totalAtoms = *totalAtoms + 14;
		break;
	default:
		throw "peptideBuilder.getGeo() wrong";
	}
	return geo;
}

vector<Geo*> initGeoChain(vector<Residue>& residuesData){
	vector<Geo*> geosData;
	int totalAtoms = 1;
	int length = residuesData.size();
	for(int i=0; i<length; i++){
		Geo* geo = getGeo(residuesData[i].m_resid,residuesData[i].m_resname, &totalAtoms);
		geo->copyMainChain(residuesData[i]);
		geo->addCB(residuesData[i]);
		geosData.push_back(geo);
	}
	return geosData;
}

vector<Atom> rebuildSideChain(vector<Geo*>& geosData, vector<Residue>& residuesData){
	vector<Atom> atomsData;
	int length = geosData.size();
	for(int i=0; i<length; i++){
		//cout << residuesData[i].m_resid << "\t" << residuesData[i].m_resname << endl;
		atomsData.push_back(residuesData[i].getAtom("N"));
		atomsData.push_back(residuesData[i].getAtom("CA"));
		atomsData.push_back(residuesData[i].getAtom("C"));
		atomsData.push_back(residuesData[i].getAtom("O"));
		if(residuesData[i].m_resname != 'G'){
			atomsData.push_back(residuesData[i].getAtom("CB"));
			if(residuesData[i].m_resname != 'A'){
				int choose_rl = residuesData[i].m_rl_id;
				int atom_length = residuesData[i].m_rotamerlibs[choose_rl].m_atoms.size();
				for(int j=0; j < atom_length; j++){
					atomsData.push_back(residuesData[i].m_rotamerlibs[choose_rl].m_atoms[j]);
				}
			}
		}

	}
	return atomsData;
}

vector<Atom> rebuildMainChainCb(vector<Geo*>& geosData, vector<Residue>& residuesData){
	vector<Atom> atomsData;
	int length = geosData.size();
	for (int i = 0; i<length; i++){
		//cout << residuesData[i].m_resid << "\t" << residuesData[i].m_resname << endl;
		atomsData.push_back(residuesData[i].getAtom("N"));
		atomsData.push_back(residuesData[i].getAtom("CA"));
		atomsData.push_back(residuesData[i].getAtom("C"));
		atomsData.push_back(residuesData[i].getAtom("O"));
		if (residuesData[i].m_resname != 'G'){
			atomsData.push_back(residuesData[i].getAtom("CB"));
		}
	}
	return atomsData;
}

vector<Atom> outputAtomsData(const vector<Geo*>& geosData){
	vector<Atom> atomsData;
	int length = geosData.size();
	for(int i=0; i<length; i++){
		geosData[i]->output(atomsData);
	}
	return atomsData;
}





/*
int main(){
	Geo* geo = getGeo('E');
	return 0;
}*/