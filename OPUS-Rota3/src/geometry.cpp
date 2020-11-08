#include <iostream> 
#include <vector> 
#include "atom.h"
#include "angle.h"
#include "residue.h"
#include "geometry.h"

using namespace std;

void Geo::addCB(Residue& residue){
}

void Geo::addSideChain(){
}

void Geo::storeSideChain(RotamerLib& rotamerlib){
}

void Geo::setRotamers(const vector<double>& rotamers){
}

void Geo::copyMainChain(Residue& residue){
	N = Atom(m_atomid, "N", m_resname, m_resid, residue.getAtom("N").m_position);
	CA = Atom(m_atomid+1, "CA", m_resname, m_resid, residue.getAtom("CA").m_position);
	C = Atom(m_atomid+2, "C", m_resname, m_resid, residue.getAtom("C").m_position);
	O = Atom(m_atomid+3, "O", m_resname, m_resid, residue.getAtom("O").m_position);
}

void Geo::output(vector<Atom>& atomsData){
	atomsData.push_back(N);
	atomsData.push_back(CA);
	atomsData.push_back(C);
	atomsData.push_back(O);
}

void GeoContainCB::addCB(Residue& residue){
	Vector3d cb = calCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle);
	CB = Atom(m_atomid+4, "CB", m_resname, m_resid, cb);
	residue.m_atoms["CB"] = CB;
}

void GeoContainCB::output(vector<Atom>& atomsData){
	Geo::output(atomsData);
	atomsData.push_back(CB);
}

GlyGeo::GlyGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 110.8914;

	C_O_length = 1.23;
	CA_C_O_angle = 120.5117;
	N_CA_C_O_diangle = 180.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	peptide_bond = 1.33;
	CA_C_N_angle = 116.642992978143;
	C_N_CA_angle = 121.382215820277;

	m_resname = 'G';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 2.5;
	m_num_dasfs = 0;
}

AlaGeo::AlaGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 111.068;

	C_O_length = 1.23;
	CA_C_O_angle = 120.5;
	N_CA_C_O_diangle = -60.5;

	phi = -120;
	psi = 140;
	omega = 180.0;
	peptide_bond = 1.33;
	CA_C_N_angle = 116.642992978143;
	C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.6860;

	m_resname = 'A';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 3.7;
	m_num_dasfs = 0;
}

SerGeo::SerGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 111.2812;

	C_O_length = 1.23;
	CA_C_O_angle = 120.5;
	N_CA_C_O_diangle = -60.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	peptide_bond = 1.33;
	CA_C_N_angle = 116.642992978143;
	C_N_CA_angle = 121.382215820277;
    

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.6618;

	CB_OG_length = 1.417;
	CA_CB_OG_angle = 110.773;
	N_CA_CB_OG_diangle = -63.3;

	m_resname = 'S';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 3.7;
	m_num_dasfs = 1;
}
void SerGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_OG_diangle = rotamers[0];
	}catch(...){
		cout << "SerGeo rotamers list: not long enough" << endl;
	}
}
void SerGeo::addSideChain(){
	Vector3d og = calCoordinates(N, CA, CB, CB_OG_length, CA_CB_OG_angle, N_CA_CB_OG_diangle);
	OG = Atom(m_atomid+5, "OG", m_resname, m_resid, og);
}
void SerGeo::storeSideChain(RotamerLib& rotamerlib){
	rotamerlib.m_atoms.push_back(OG);
}
void SerGeo::output(vector<Atom>& atomsData){
	GeoContainCB::output(atomsData);
	atomsData.push_back(OG);
}

CysGeo::CysGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 110.8856;

	C_O_length = 1.23;
	CA_C_O_angle = 120.5;
	N_CA_C_O_diangle = -60.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	peptide_bond = 1.33;
	CA_C_N_angle = 116.642992978143;
	C_N_CA_angle = 121.382215820277;


	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.5037;

	CB_SG_length = 1.808;
	CA_CB_SG_angle = 113.8169;
	N_CA_CB_SG_diangle = -62.2;

	m_resname = 'C';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 3.4;
	m_num_dasfs = 1;
}
void CysGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_SG_diangle = rotamers[0];
	}catch(...){
		cout << "Geo rotamers list: not long enough" << endl;
	}
}
void CysGeo::addSideChain(){
	Vector3d sg = calCoordinates(N, CA, CB, CB_SG_length, CA_CB_SG_angle, N_CA_CB_SG_diangle);
	SG = Atom(m_atomid+5, "SG", m_resname, m_resid, sg);
}
void CysGeo::storeSideChain(RotamerLib& rotamerlib){
	rotamerlib.m_atoms.push_back(SG);
}
void CysGeo::output(vector<Atom>& atomsData){
	GeoContainCB::output(atomsData);
	atomsData.push_back(SG);
}
                                       
ValGeo::ValGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 109.7698;
        
	C_O_length = 1.23;
	CA_C_O_angle = 120.5686;
	N_CA_C_O_diangle = -60.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	peptide_bond = 1.33;
	CA_C_N_angle = 116.642992978143;
	C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 123.2347;

	CB_CG1_length = 1.527;
	CA_CB_CG1_angle = 110.7;
	N_CA_CB_CG1_diangle = 177.2;

	CB_CG2_length = 1.527;
	CG1_CB_CG2_angle = 109.50; 
	CA_CG1_CB_CG2_diangle = -120;

	m_resname = 'V';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 3.5;
	m_num_dasfs = 1;
}
void ValGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG1_diangle = rotamers[0];
	}catch(...){
		cout << "ValGeo rotamers list: not long enough" << endl;
	}
}
void ValGeo::addSideChain(){
	Vector3d cg1 = calCoordinates(N, CA, CB, CB_CG1_length, CA_CB_CG1_angle, N_CA_CB_CG1_diangle);
	CG1 = Atom(m_atomid+5, "CG1", m_resname, m_resid, cg1);
	Vector3d cg2 = calCoordinates(CA, CG1, CB, CB_CG2_length, CG1_CB_CG2_angle, CA_CG1_CB_CG2_diangle);
	CG2 = Atom(m_atomid+6, "CG2", m_resname, m_resid, cg2);
}
void ValGeo::storeSideChain(RotamerLib& rotamerlib){
	rotamerlib.m_atoms.push_back(CG1);
	rotamerlib.m_atoms.push_back(CG2);
}
void ValGeo::output(vector<Atom>& atomsData){
	GeoContainCB::output(atomsData);
	atomsData.push_back(CG1);
	atomsData.push_back(CG2);
}

IleGeo::IleGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 109.7202;

	C_O_length = 1.23;
	CA_C_O_angle = 120.5403;
	N_CA_C_O_diangle = -60.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	peptide_bond = 1.33;
	CA_C_N_angle = 116.642992978143;
	C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 123.2347;

	CB_CG1_length = 1.527;
	CA_CB_CG1_angle = 110.7;
	N_CA_CB_CG1_diangle = 59.7;

	CB_CG2_length = 1.527;
	CG1_CB_CG2_angle = 109.50;
	CA_CG1_CB_CG2_diangle = 120;

	CG1_CD1_length = 1.52;
	CB_CG1_CD1_angle = 113.97;
	CA_CB_CG1_CD1_diangle = 169.8;

	m_resname = 'I';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 3.5;
	m_num_dasfs = 2;
}
void IleGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG1_diangle = rotamers[0];
		CA_CB_CG1_CD1_diangle = rotamers[1];
	}catch(...){
		cout << "IleGeo rotamers list: not long enough" << endl;
	}
}
void IleGeo::addSideChain(){
	Vector3d cg1 = calCoordinates(N, CA, CB, CB_CG1_length, CA_CB_CG1_angle, N_CA_CB_CG1_diangle);
	CG1 = Atom(m_atomid+5, "CG1", m_resname, m_resid, cg1);
	Vector3d cg2 = calCoordinates(CA, CG1, CB, CB_CG2_length, CG1_CB_CG2_angle, CA_CG1_CB_CG2_diangle);
	CG2 = Atom(m_atomid+6, "CG2", m_resname, m_resid, cg2);
	Vector3d cd1 = calCoordinates(CA, CB, CG1, CG1_CD1_length, CB_CG1_CD1_angle, CA_CB_CG1_CD1_diangle);
	CD1 = Atom(m_atomid+7, "CD1", m_resname, m_resid, cd1);
}
void IleGeo::storeSideChain(RotamerLib& rotamerlib){
	rotamerlib.m_atoms.push_back(CG1);
	rotamerlib.m_atoms.push_back(CG2);
	rotamerlib.m_atoms.push_back(CD1);
}
void IleGeo::output(vector<Atom>& atomsData){
	GeoContainCB::output(atomsData);
	atomsData.push_back(CG1);
	atomsData.push_back(CG2);
	atomsData.push_back(CD1);
}

LeuGeo::LeuGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 110.8652;

	C_O_length = 1.23;
	CA_C_O_angle = 120.4647;
	N_CA_C_O_diangle = 120.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	peptide_bond = 1.33;
	CA_C_N_angle = 116.642992978143;
	C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.4948;

	CB_CG_length = 1.53;
	CA_CB_CG_angle = 116.10;
	N_CA_CB_CG_diangle = -60.1;

	CG_CD1_length = 1.524;
	CB_CG_CD1_angle = 110.27;
	CA_CB_CG_CD1_diangle = 174.9;

	CG_CD2_length = 1.525;
	CD1_CG_CD2_angle = 109.50;
	CB_CD1_CG_CD2_diangle = -120;

	m_resname =  'L';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 3.6;
	m_num_dasfs = 2;
}
void LeuGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG_diangle = rotamers[0];
		CA_CB_CG_CD1_diangle = rotamers[1];
	}catch(...){
		cout << "LeuGeo rotamers list: not long enough" << endl;
	}
}
void LeuGeo::addSideChain(){
	Vector3d cg = calCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
	CG = Atom(m_atomid+5, "CG", m_resname, m_resid, cg);
	Vector3d cd1 = calCoordinates(CA, CB, CG, CG_CD1_length, CB_CG_CD1_angle, CA_CB_CG_CD1_diangle);
	CD1 = Atom(m_atomid+6, "CD1", m_resname, m_resid, cd1);
	Vector3d cd2 = calCoordinates(CB, CD1, CG, CG_CD2_length, CD1_CG_CD2_angle, CB_CD1_CG_CD2_diangle);
	CD2 = Atom(m_atomid+7, "CD2", m_resname, m_resid, cd2);
}
void LeuGeo::storeSideChain(RotamerLib& rotamerlib){
	rotamerlib.m_atoms.push_back(CG);
	rotamerlib.m_atoms.push_back(CD1);
	rotamerlib.m_atoms.push_back(CD2);
}
void LeuGeo::output(vector<Atom>& atomsData){
	GeoContainCB::output(atomsData);
	atomsData.push_back(CG);
	atomsData.push_back(CD1);
	atomsData.push_back(CD2);
}

ThrGeo::ThrGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 110.7014;

	C_O_length = 1.23;
	CA_C_O_angle = 120.5359;
	N_CA_C_O_diangle = 120.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	peptide_bond = 1.33;
	CA_C_N_angle = 116.642992978143;
	C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 123.0953;

	CB_OG1_length = 1.43;
	CA_CB_OG1_angle = 109.18;
	N_CA_CB_OG1_diangle = 60.0;

	CB_CG2_length = 1.53;
	OG1_CB_CG2_angle = 109.50;
	CA_OG1_CB_CG2_diangle = 120;

	m_resname =  'T';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 3.6;
	m_num_dasfs = 1;
}
void ThrGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_OG1_diangle = rotamers[0];
	}catch(...){
		cout << "ThrGeo rotamers list: not long enough" << endl;
	}
}
void ThrGeo::addSideChain(){
	Vector3d og1 = calCoordinates(N, CA, CB, CB_OG1_length, CA_CB_OG1_angle, N_CA_CB_OG1_diangle);
	OG1 = Atom(m_atomid+5, "OG1", m_resname, m_resid, og1);
	Vector3d cg2 = calCoordinates(CA, OG1, CB, CB_CG2_length, OG1_CB_CG2_angle, CA_OG1_CB_CG2_diangle);
	CG2 = Atom(m_atomid+6, "CG2", m_resname, m_resid, cg2);
}
void ThrGeo::storeSideChain(RotamerLib& rotamerlib){
	rotamerlib.m_atoms.push_back(OG1);
	rotamerlib.m_atoms.push_back(CG2);
}
void ThrGeo::output(vector<Atom>& atomsData){
	GeoContainCB::output(atomsData);
	atomsData.push_back(OG1);
	atomsData.push_back(CG2);
}

ArgGeo::ArgGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 110.98;

	C_O_length = 1.23;
	CA_C_O_angle = 120.54;
	N_CA_C_O_diangle = 120.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	peptide_bond = 1.33;
	CA_C_N_angle = 116.642992978143;
	C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.76;

	CB_CG_length = 1.52;
	CA_CB_CG_angle = 113.83;
	N_CA_CB_CG_diangle = -65.2;

	CG_CD_length = 1.52;
	CB_CG_CD_angle = 111.79;
	CA_CB_CG_CD_diangle = -179.2;

	CD_NE_length = 1.46;
	CG_CD_NE_angle = 111.68;
	CB_CG_CD_NE_diangle = -179.3;

	NE_CZ_length = 1.33;
	CD_NE_CZ_angle = 124.79;
	CG_CD_NE_CZ_diangle = -178.7;

	CZ_NH1_length = 1.33;
	NE_CZ_NH1_angle = 120.64;
	CD_NE_CZ_NH1_diangle = 0.0;

	CZ_NH2_length = 1.33;
	NE_CZ_NH2_angle = 119.63;
	CD_NE_CZ_NH2_diangle = 180.0;

	m_resname =  'R';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 6.2;
	m_num_dasfs = 4;
}
void ArgGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG_diangle = rotamers[0];
		CA_CB_CG_CD_diangle = rotamers[1];
		CB_CG_CD_NE_diangle = rotamers[2];
		CG_CD_NE_CZ_diangle = rotamers[3];
	}catch(...){
		cout << "ArgGeo rotamers list: not long enough" << endl;
	}
}
void ArgGeo::addSideChain(){
	Vector3d cg = calCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
	CG = Atom(m_atomid+5, "CG", m_resname, m_resid, cg);
	Vector3d cd = calCoordinates(CA, CB, CG, CG_CD_length, CB_CG_CD_angle, CA_CB_CG_CD_diangle);
	CD = Atom(m_atomid+6, "CD", m_resname, m_resid, cd);
	Vector3d ne = calCoordinates(CB, CG, CD, CD_NE_length, CG_CD_NE_angle, CB_CG_CD_NE_diangle);
	NE = Atom(m_atomid+7, "NE", m_resname, m_resid, ne);
	Vector3d cz = calCoordinates(CG, CD, NE, NE_CZ_length, CD_NE_CZ_angle, CG_CD_NE_CZ_diangle);
	CZ = Atom(m_atomid+8, "CZ", m_resname, m_resid, cz);
	Vector3d nh1 = calCoordinates(CD, NE, CZ, CZ_NH1_length, NE_CZ_NH1_angle, CD_NE_CZ_NH1_diangle);
	NH1 = Atom(m_atomid+9, "NH1", m_resname, m_resid, nh1);
	Vector3d nh2 = calCoordinates(CD, NE, CZ, CZ_NH2_length, NE_CZ_NH2_angle, CD_NE_CZ_NH2_diangle);
	NH2 = Atom(m_atomid+10, "NH2", m_resname, m_resid, nh2);
}
void ArgGeo::storeSideChain(RotamerLib& rotamerlib){
	rotamerlib.m_atoms.push_back(CG);
	rotamerlib.m_atoms.push_back(CD);
	rotamerlib.m_atoms.push_back(NE);
	rotamerlib.m_atoms.push_back(CZ);
	rotamerlib.m_atoms.push_back(NH1);
	rotamerlib.m_atoms.push_back(NH2);
}
void ArgGeo::output(vector<Atom>& atomsData){
	GeoContainCB::output(atomsData);
	atomsData.push_back(CG);
	atomsData.push_back(CD);
	atomsData.push_back(NE);
	atomsData.push_back(CZ);
	atomsData.push_back(NH1);
	atomsData.push_back(NH2);
}

LysGeo::LysGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 111.08;

	C_O_length = 1.23;
	CA_C_O_angle = 120.54;
	N_CA_C_O_diangle = 120.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	peptide_bond = 1.33;
	CA_C_N_angle = 116.642992978143;
	C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.76;

	CB_CG_length = 1.52;
	CA_CB_CG_angle = 113.83;
	N_CA_CB_CG_diangle = -64.5;

	CG_CD_length = 1.52;
	CB_CG_CD_angle = 111.79;
	CA_CB_CG_CD_diangle = -178.1;

	CD_CE_length = 1.46;
	CG_CD_CE_angle = 111.68;
	CB_CG_CD_CE_diangle = -179.6;

	CE_NZ_length = 1.33;
	CD_CE_NZ_angle = 124.79;
	CG_CD_CE_NZ_diangle = 179.6;

	m_resname =  'K';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 5;
	m_num_dasfs = 4;
}
void LysGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG_diangle = rotamers[0];
		CA_CB_CG_CD_diangle = rotamers[1];
		CB_CG_CD_CE_diangle = rotamers[2];
		CG_CD_CE_NZ_diangle = rotamers[3];
	}catch(...){
		cout << "LysGeo rotamers list: not long enough" << endl;
	}
}
void LysGeo::addSideChain(){
	Vector3d cg = calCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
	CG = Atom(m_atomid+5, "CG", m_resname, m_resid, cg);
	Vector3d cd = calCoordinates(CA, CB, CG, CG_CD_length, CB_CG_CD_angle, CA_CB_CG_CD_diangle);
	CD = Atom(m_atomid+6, "CD", m_resname, m_resid, cd);
	Vector3d ce = calCoordinates(CB, CG, CD, CD_CE_length, CG_CD_CE_angle, CB_CG_CD_CE_diangle);
	CE = Atom(m_atomid+7, "CE", m_resname, m_resid, ce);
	Vector3d nz = calCoordinates(CG, CD, CE, CE_NZ_length, CD_CE_NZ_angle, CG_CD_CE_NZ_diangle);
	NZ = Atom(m_atomid+8, "NZ", m_resname, m_resid, nz);
}
void LysGeo::storeSideChain(RotamerLib& rotamerlib){
	rotamerlib.m_atoms.push_back(CG);
	rotamerlib.m_atoms.push_back(CD);
	rotamerlib.m_atoms.push_back(CE);
	rotamerlib.m_atoms.push_back(NZ);
}
void LysGeo::output(vector<Atom>& atomsData){
	GeoContainCB::output(atomsData);
	atomsData.push_back(CG);
	atomsData.push_back(CD);
	atomsData.push_back(CE);
	atomsData.push_back(NZ);

}

AspGeo::AspGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 111.03;    

	C_O_length = 1.23;
	CA_C_O_angle = 120.51;
	N_CA_C_O_diangle = 120.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	peptide_bond = 1.33;
	CA_C_N_angle = 116.642992978143;
	C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.82;

	CB_CG_length = 1.52;
	CA_CB_CG_angle = 113.06;
	N_CA_CB_CG_diangle = -66.4;

	CG_OD1_length = 1.25;
	CB_CG_OD1_angle = 119.22;
	CA_CB_CG_OD1_diangle = -46.7;

	CG_OD2_length = 1.25;
	CB_CG_OD2_angle = 118.218;
	CA_CB_CG_OD2_diangle = 180+CA_CB_CG_OD1_diangle;

	m_resname =  'D';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 3.6;
	m_num_dasfs = 2;
}
void AspGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG_diangle = rotamers[0];
		CA_CB_CG_OD1_diangle = rotamers[1];
		if (CA_CB_CG_OD1_diangle > 0){
			CA_CB_CG_OD2_diangle = CA_CB_CG_OD1_diangle-180.0;
		}else{
			CA_CB_CG_OD2_diangle = CA_CB_CG_OD1_diangle+180.0;
		}
	}catch(...){
		cout << "AspGeo rotamers list: not long enough" << endl;
	}
}
void AspGeo::addSideChain(){
	Vector3d cg = calCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
	CG = Atom(m_atomid+5, "CG", m_resname, m_resid, cg);
	Vector3d od1 = calCoordinates(CA, CB, CG, CG_OD1_length, CB_CG_OD1_angle, CA_CB_CG_OD1_diangle);
	OD1 = Atom(m_atomid+6, "OD1", m_resname, m_resid, od1);
	Vector3d od2 = calCoordinates(CA, CB, CG, CG_OD2_length, CB_CG_OD2_angle, CA_CB_CG_OD2_diangle);
	OD2 = Atom(m_atomid+7, "OD2", m_resname, m_resid, od2);
}
void AspGeo::storeSideChain(RotamerLib& rotamerlib){
	rotamerlib.m_atoms.push_back(CG);
	rotamerlib.m_atoms.push_back(OD1);
	rotamerlib.m_atoms.push_back(OD2);
}
void AspGeo::output(vector<Atom>& atomsData){
	GeoContainCB::output(atomsData);
	atomsData.push_back(CG);
	atomsData.push_back(OD1);
	atomsData.push_back(OD2);
}

AsnGeo::AsnGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 111.5;

	C_O_length = 1.23;
	CA_C_O_angle = 120.4826;
	N_CA_C_O_diangle =  -60.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	peptide_bond = 1.33;
	CA_C_N_angle = 116.642992978143;
	C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 123.2254;

	CB_CG_length = 1.52;
	CA_CB_CG_angle = 112.62;
	N_CA_CB_CG_diangle = -65.5;

	CG_OD1_length = 1.23;
	CB_CG_OD1_angle = 120.85;
	CA_CB_CG_OD1_diangle = -58.3;

	CG_ND2_length = 1.33;
	CB_CG_ND2_angle = 116.48;
	CA_CB_CG_ND2_diangle = 180.0+CA_CB_CG_OD1_diangle;

	m_resname =  'N';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 3.6;
	m_num_dasfs = 2;
}
void AsnGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG_diangle = rotamers[0];
		CA_CB_CG_OD1_diangle = rotamers[1];
		if (CA_CB_CG_OD1_diangle > 0){
			CA_CB_CG_ND2_diangle = CA_CB_CG_OD1_diangle-180.0;
		}else{
			CA_CB_CG_ND2_diangle = CA_CB_CG_OD1_diangle+180.0;
		}
	}catch(...){
		cout << "AsnGeo rotamers list: not long enough" << endl;
	}
}
void AsnGeo::addSideChain(){
	Vector3d cg = calCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
	CG = Atom(m_atomid+5, "CG", m_resname, m_resid, cg);
	Vector3d od1 = calCoordinates(CA, CB, CG, CG_OD1_length, CB_CG_OD1_angle, CA_CB_CG_OD1_diangle);
	OD1 = Atom(m_atomid+6, "OD1", m_resname, m_resid, od1);
	Vector3d nd2 = calCoordinates(CA, CB, CG, CG_ND2_length, CB_CG_ND2_angle, CA_CB_CG_ND2_diangle);
	ND2 = Atom(m_atomid+7, "ND2", m_resname, m_resid, nd2);
}
void AsnGeo::storeSideChain(RotamerLib& rotamerlib){
	rotamerlib.m_atoms.push_back(CG);
	rotamerlib.m_atoms.push_back(OD1);
	rotamerlib.m_atoms.push_back(ND2);
}
void AsnGeo::output(vector<Atom>& atomsData){
	GeoContainCB::output(atomsData);
	atomsData.push_back(CG);
	atomsData.push_back(OD1);
	atomsData.push_back(ND2);
}

GluGeo::GluGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 111.1703;

	C_O_length = 1.23;
	CA_C_O_angle = 120.511;
	N_CA_C_O_diangle = 120.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	peptide_bond = 1.33;
	CA_C_N_angle = 116.642992978143;
	C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.8702;

	CB_CG_length = 1.52;
	CA_CB_CG_angle = 113.82;
	N_CA_CB_CG_diangle = -63.8;

	CG_CD_length = 1.52;
	CB_CG_CD_angle = 113.31;
	CA_CB_CG_CD_diangle = -179.8;

	CD_OE1_length = 1.25;
	CG_CD_OE1_angle = 119.02;
	CB_CG_CD_OE1_diangle = -6.2;

	CD_OE2_length = 1.25;
	CG_CD_OE2_angle = 118.08;
	CB_CG_CD_OE2_diangle = 180.0+CB_CG_CD_OE1_diangle;

	m_resname =  'E';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 3.7;
	m_num_dasfs = 3;
}
void GluGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG_diangle = rotamers[0];
		CA_CB_CG_CD_diangle = rotamers[1];
		CB_CG_CD_OE1_diangle = rotamers[2];
		if (CB_CG_CD_OE1_diangle > 0){
			CB_CG_CD_OE2_diangle = CB_CG_CD_OE1_diangle-180.0;
		}else{
			CB_CG_CD_OE2_diangle = CB_CG_CD_OE1_diangle+180.0;
		}
	}catch(...){
		cout << "GluGeo rotamers list: not long enough" << endl;
	}
}
void GluGeo::addSideChain(){
	Vector3d cg = calCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
	CG = Atom(m_atomid+5, "CG", m_resname, m_resid, cg);
	Vector3d cd = calCoordinates(CA, CB, CG, CG_CD_length, CB_CG_CD_angle, CA_CB_CG_CD_diangle);
	CD = Atom(m_atomid+6, "CD", m_resname, m_resid, cd);
	Vector3d oe1 = calCoordinates(CB, CG, CD, CD_OE1_length, CG_CD_OE1_angle, CB_CG_CD_OE1_diangle);
	OE1 = Atom(m_atomid+7, "OE1", m_resname, m_resid, oe1);
	Vector3d oe2 = calCoordinates(CB, CG, CD, CD_OE2_length, CG_CD_OE2_angle, CB_CG_CD_OE2_diangle);
	OE2 = Atom(m_atomid+8, "OE2", m_resname, m_resid, oe2);
}
void GluGeo::storeSideChain(RotamerLib& rotamerlib){
	rotamerlib.m_atoms.push_back(CG);
	rotamerlib.m_atoms.push_back(CD);
	rotamerlib.m_atoms.push_back(OE1);
	rotamerlib.m_atoms.push_back(OE2);
}
void GluGeo::output(vector<Atom>& atomsData){
	GeoContainCB::output(atomsData);
	atomsData.push_back(CG);
	atomsData.push_back(CD);
	atomsData.push_back(OE1);
	atomsData.push_back(OE2);
}

GlnGeo::GlnGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 111.0849;

	C_O_length = 1.23;
	CA_C_O_angle = 120.5029;
	N_CA_C_O_diangle = 120.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	peptide_bond = 1.33;
	CA_C_N_angle = 116.642992978143;
	C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.8134;

	CB_CG_length = 1.52;
	CA_CB_CG_angle = 113.75;
	N_CA_CB_CG_diangle = -60.2;

	CG_CD_length = 1.52;
	CB_CG_CD_angle = 112.78;
	CA_CB_CG_CD_diangle = -69.6;

	CD_OE1_length = 1.24;
	CG_CD_OE1_angle = 120.86;
	CB_CG_CD_OE1_diangle = -50.5;

	CD_NE2_length = 1.33;
	CG_CD_NE2_angle = 116.50;
	CB_CG_CD_NE2_diangle = 180+CB_CG_CD_OE1_diangle;

	m_resname =  'Q';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 3.8;
	m_num_dasfs = 3;
}
void GlnGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG_diangle = rotamers[0];
		CA_CB_CG_CD_diangle = rotamers[1];
		CB_CG_CD_OE1_diangle = rotamers[2];
		if (CB_CG_CD_OE1_diangle > 0){
			CB_CG_CD_NE2_diangle = CB_CG_CD_OE1_diangle-180.0;
		}else{
			CB_CG_CD_NE2_diangle = CB_CG_CD_OE1_diangle+180.0;
		}
	}catch(...){
		cout << "GlnGeo rotamers list: not long enough" << endl;
	}
}
void GlnGeo::addSideChain(){
	Vector3d cg = calCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
	CG = Atom(m_atomid+5, "CG", m_resname, m_resid, cg);
	Vector3d cd = calCoordinates(CA, CB, CG, CG_CD_length, CB_CG_CD_angle, CA_CB_CG_CD_diangle);
	CD = Atom(m_atomid+6, "CD", m_resname, m_resid, cd);
	Vector3d oe1 = calCoordinates(CB, CG, CD, CD_OE1_length, CG_CD_OE1_angle, CB_CG_CD_OE1_diangle);
	OE1 = Atom(m_atomid+7, "OE1", m_resname, m_resid, oe1);
	Vector3d ne2 = calCoordinates(CB, CG, CD, CD_NE2_length, CG_CD_NE2_angle, CB_CG_CD_NE2_diangle);
	NE2 = Atom(m_atomid+8, "NE2", m_resname, m_resid, ne2);
}
void GlnGeo::storeSideChain(RotamerLib& rotamerlib){
	rotamerlib.m_atoms.push_back(CG);
	rotamerlib.m_atoms.push_back(CD);
	rotamerlib.m_atoms.push_back(OE1);
	rotamerlib.m_atoms.push_back(NE2);
}
void GlnGeo::output(vector<Atom>& atomsData){
	GeoContainCB::output(atomsData);
	atomsData.push_back(CG);
	atomsData.push_back(CD);
	atomsData.push_back(OE1);
	atomsData.push_back(NE2);
}

MetGeo::MetGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 110.9416;

	C_O_length = 1.23;
	CA_C_O_angle = 120.4816;
	N_CA_C_O_diangle = 120.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	peptide_bond = 1.33;
	CA_C_N_angle = 116.642992978143;
	C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.6733;

	CB_CG_length = 1.52;
	CA_CB_CG_angle =  113.68;
	N_CA_CB_CG_diangle = -64.4;

	CG_SD_length = 1.81;
	CB_CG_SD_angle = 112.69;
	CA_CB_CG_SD_diangle = -179.6;

	SD_CE_length = 1.79;
	CG_SD_CE_angle = 100.61;
	CB_CG_SD_CE_diangle = 70.1;

	m_resname = 'M';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 4.2;
	m_num_dasfs = 3;

}
void MetGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG_diangle = rotamers[0];
		CA_CB_CG_SD_diangle = rotamers[1];
		CB_CG_SD_CE_diangle = rotamers[2];
	}catch(...){
		cout << "MetGeo rotamers list: not long enough" << endl;
	}
}
void MetGeo::addSideChain(){
	Vector3d cg = calCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
	CG = Atom(m_atomid+5, "CG", m_resname, m_resid, cg);
	Vector3d sd = calCoordinates(CA, CB, CG, CG_SD_length, CB_CG_SD_angle, CA_CB_CG_SD_diangle);
	SD = Atom(m_atomid+6, "SD", m_resname, m_resid, sd);
	Vector3d ce = calCoordinates(CB, CG, SD, SD_CE_length, CG_SD_CE_angle, CB_CG_SD_CE_diangle);
	CE = Atom(m_atomid+7, "CE", m_resname, m_resid, ce);
}
void MetGeo::storeSideChain(RotamerLib& rotamerlib){
	rotamerlib.m_atoms.push_back(CG);
	rotamerlib.m_atoms.push_back(SD);
	rotamerlib.m_atoms.push_back(CE);
}
void MetGeo::output(vector<Atom>& atomsData){
	GeoContainCB::output(atomsData);
	atomsData.push_back(CG);
	atomsData.push_back(SD);
	atomsData.push_back(CE);
}

HisGeo::HisGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 111.0859;

	C_O_length = 1.23;
	CA_C_O_angle = 120.4732;
	N_CA_C_O_diangle = 120.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	peptide_bond = 1.33;
	CA_C_N_angle = 116.642992978143;
	C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.6711;

	CB_CG_length = 1.49;
	CA_CB_CG_angle = 113.74;
	N_CA_CB_CG_diangle = -63.2;

	CG_ND1_length = 1.38;
	CB_CG_ND1_angle = 122.85;
	CA_CB_CG_ND1_diangle = -75.7;    

	CG_CD2_length = 1.35;
	CB_CG_CD2_angle = 130.61;
	CA_CB_CG_CD2_diangle = 180.0+CA_CB_CG_ND1_diangle;

	ND1_CE1_length = 1.32;
	CG_ND1_CE1_angle = 108.5;
	CB_CG_ND1_CE1_diangle = 180.0;

	CD2_NE2_length = 1.35;
	CG_CD2_NE2_angle = 108.5;
	CB_CG_CD2_NE2_diangle = 180.0;

	m_resname =  'H';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 3.7;
	m_num_dasfs = 2;
}
void HisGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG_diangle = rotamers[0];
		CA_CB_CG_ND1_diangle = rotamers[1];
		if (CA_CB_CG_ND1_diangle > 0){
			CA_CB_CG_CD2_diangle = CA_CB_CG_ND1_diangle-180.0;
		}else{
			CA_CB_CG_CD2_diangle = CA_CB_CG_ND1_diangle+180.0;
		}
	}catch(...){
		cout << "HisGeo rotamers list: not long enough" << endl;
	}
}
void HisGeo::addSideChain(){
	Vector3d cg = calCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
	CG = Atom(m_atomid+5, "CG", m_resname, m_resid, cg);
	Vector3d nd1 = calCoordinates(CA, CB, CG, CG_ND1_length, CB_CG_ND1_angle, CA_CB_CG_ND1_diangle);
	ND1 = Atom(m_atomid+6, "ND1", m_resname, m_resid, nd1);
	Vector3d cd2 = calCoordinates(CA, CB, CG, CG_CD2_length, CB_CG_CD2_angle, CA_CB_CG_CD2_diangle);
	CD2 = Atom(m_atomid+7, "CD2", m_resname, m_resid, cd2);
	Vector3d ce1 = calCoordinates(CB, CG, ND1, ND1_CE1_length, CG_ND1_CE1_angle, CB_CG_ND1_CE1_diangle);
	CE1 = Atom(m_atomid+8, "CE1", m_resname, m_resid, ce1);
	Vector3d ne2 = calCoordinates(CB, CG, CD2, CD2_NE2_length, CG_CD2_NE2_angle, CB_CG_CD2_NE2_diangle);
	NE2 = Atom(m_atomid+9, "NE2", m_resname, m_resid, ne2);
}
void HisGeo::storeSideChain(RotamerLib& rotamerlib){
	rotamerlib.m_atoms.push_back(CG);
	rotamerlib.m_atoms.push_back(ND1);
	rotamerlib.m_atoms.push_back(CD2);
	rotamerlib.m_atoms.push_back(CE1);
	rotamerlib.m_atoms.push_back(NE2);
}
void HisGeo::output(vector<Atom>& atomsData){
	GeoContainCB::output(atomsData);
	atomsData.push_back(CG);
	atomsData.push_back(ND1);
	atomsData.push_back(CD2);
	atomsData.push_back(CE1);
	atomsData.push_back(NE2);
}

ProGeo::ProGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 112.7499;

	C_O_length = 1.23;
	CA_C_O_angle = 120.2945;
	N_CA_C_O_diangle = -45.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	peptide_bond = 1.33;
	CA_C_N_angle  = 116.642992978143;
	C_N_CA_angle =  121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 115.2975;

	CB_CG_length = 1.49;
	CA_CB_CG_angle = 104.21;
	N_CA_CB_CG_diangle = 29.6;

	CG_CD_length = 1.50;
	CB_CG_CD_angle = 105.03;
	CA_CB_CG_CD_diangle = -34.8;

	m_resname =  'P';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 3.5;
	m_num_dasfs = 2;
}
void ProGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG_diangle = rotamers[0];
		CA_CB_CG_CD_diangle = rotamers[1];
	}catch(...){
		cout << "ProGeo rotamers list: not long enough" << endl;
	}
}
void ProGeo::addSideChain(){
	Vector3d cg = calCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
	CG = Atom(m_atomid+5, "CG", m_resname, m_resid, cg);
	Vector3d cd = calCoordinates(CA, CB, CG, CG_CD_length, CB_CG_CD_angle, CA_CB_CG_CD_diangle);
	CD = Atom(m_atomid+6, "CD", m_resname, m_resid, cd);
}
void ProGeo::storeSideChain(RotamerLib& rotamerlib){
	rotamerlib.m_atoms.push_back(CG);
	rotamerlib.m_atoms.push_back(CD);
}
void ProGeo::output(vector<Atom>& atomsData){
	GeoContainCB::output(atomsData);
	atomsData.push_back(CG);
	atomsData.push_back(CD);
}

PheGeo::PheGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 110.7528;

	C_O_length = 1.23;
	CA_C_O_angle = 120.5316;
	N_CA_C_O_diangle = 120.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	peptide_bond = 1.33;
	CA_C_N_angle = 116.642992978143;
	C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.6054;

	CB_CG_length = 1.50;
	CA_CB_CG_angle = 113.85;
	N_CA_CB_CG_diangle = -64.7;

	CG_CD1_length = 1.39;
	CB_CG_CD1_angle = 120.0;
	CA_CB_CG_CD1_diangle = 93.3;

	CG_CD2_length = 1.39;
	CB_CG_CD2_angle = 120.0;
	CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle-180.0;

	CD1_CE1_length = 1.39;
	CG_CD1_CE1_angle = 120.0;
	CB_CG_CD1_CE1_diangle = 180.0;

	CD2_CE2_length = 1.39;
	CG_CD2_CE2_angle = 120.0;
	CB_CG_CD2_CE2_diangle = 180.0;

	CE1_CZ_length = 1.39;
	CD1_CE1_CZ_angle = 120.0;
	CG_CD1_CE1_CZ_diangle = 0.0;

	m_resname =  'F';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 4.3;
	m_num_dasfs = 2;
}
void PheGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG_diangle = rotamers[0];
		CA_CB_CG_CD1_diangle = rotamers[1];
		if (CA_CB_CG_CD1_diangle > 0){
			CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle-180.0;
		}else{
			CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle+180.0;
		}
	}catch(...){
		cout << "PheGeo rotamers list: not long enough" << endl;
	}
}
void PheGeo::addSideChain(){
	Vector3d cg = calCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
	CG = Atom(m_atomid+5, "CG", m_resname, m_resid, cg);
	Vector3d cd1 = calCoordinates(CA, CB, CG, CG_CD1_length, CB_CG_CD1_angle, CA_CB_CG_CD1_diangle);
	CD1 = Atom(m_atomid+6, "CD1", m_resname, m_resid, cd1);
	Vector3d cd2 = calCoordinates(CA, CB, CG, CG_CD2_length, CB_CG_CD2_angle, CA_CB_CG_CD2_diangle);
	CD2 = Atom(m_atomid+7, "CD2", m_resname, m_resid, cd2);
	Vector3d ce1 = calCoordinates(CB, CG, CD1, CD1_CE1_length, CG_CD1_CE1_angle, CB_CG_CD1_CE1_diangle);
	CE1 = Atom(m_atomid+8, "CE1", m_resname, m_resid, ce1);
	Vector3d ce2 = calCoordinates(CB, CG, CD2, CD2_CE2_length, CG_CD2_CE2_angle, CB_CG_CD2_CE2_diangle);
	CE2 = Atom(m_atomid+9, "CE2", m_resname, m_resid, ce2);
	Vector3d cz = calCoordinates(CG, CD1, CE1, CE1_CZ_length, CD1_CE1_CZ_angle, CG_CD1_CE1_CZ_diangle);
	CZ = Atom(m_atomid+10, "CZ", m_resname, m_resid, cz);
}
void PheGeo::storeSideChain(RotamerLib& rotamerlib){
	rotamerlib.m_atoms.push_back(CG);
	rotamerlib.m_atoms.push_back(CD1);
	rotamerlib.m_atoms.push_back(CD2);
	rotamerlib.m_atoms.push_back(CE1);
	rotamerlib.m_atoms.push_back(CE2);
	rotamerlib.m_atoms.push_back(CZ);
}
void PheGeo::output(vector<Atom>& atomsData){
	GeoContainCB::output(atomsData);
	atomsData.push_back(CG);
	atomsData.push_back(CD1);
	atomsData.push_back(CD2);
	atomsData.push_back(CE1);
	atomsData.push_back(CE2);
	atomsData.push_back(CZ);
}

TyrGeo::TyrGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 110.9288;

	C_O_length = 1.23;
	CA_C_O_angle = 120.5434;
	N_CA_C_O_diangle = 120.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	peptide_bond = 1.33;
	CA_C_N_angle = 116.642992978143;
	C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.6023;

	CB_CG_length = 1.51;
	CA_CB_CG_angle =  113.8;
	N_CA_CB_CG_diangle = -64.3;

	CG_CD1_length = 1.39;
	CB_CG_CD1_angle = 120.98;
	CA_CB_CG_CD1_diangle = 93.1;

	CG_CD2_length = 1.39;
	CB_CG_CD2_angle = 120.82;
	CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle-180.0;

	CD1_CE1_length = 1.39;
	CG_CD1_CE1_angle = 120.0;
	CB_CG_CD1_CE1_diangle = 180.0;

	CD2_CE2_length = 1.39;
	CG_CD2_CE2_angle = 120.0;
	CB_CG_CD2_CE2_diangle = 180.0;

	CE1_CZ_length = 1.39;
	CD1_CE1_CZ_angle = 120.0;
	CG_CD1_CE1_CZ_diangle = 0.0;

	CZ_OH_length = 1.39;
	CE1_CZ_OH_angle = 119.78;
	CD1_CE1_CZ_OH_diangle = 180.0;

	m_resname =  'Y';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 5.7;
	m_num_dasfs = 2;
}
void TyrGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG_diangle = rotamers[0];
		CA_CB_CG_CD1_diangle = rotamers[1];
		if (CA_CB_CG_CD1_diangle > 0){
			CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle-180.0;
		}else{
			CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle+180.0;
		}
	}catch(...){
		cout << "TyrGeo rotamers list: not long enough" << endl;
	}

}
void TyrGeo::addSideChain(){
	Vector3d cg = calCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
	CG = Atom(m_atomid+5, "CG", m_resname, m_resid, cg);
	Vector3d cd1 = calCoordinates(CA, CB, CG, CG_CD1_length, CB_CG_CD1_angle, CA_CB_CG_CD1_diangle);
	CD1 = Atom(m_atomid+6, "CD1", m_resname, m_resid, cd1);
	Vector3d cd2 = calCoordinates(CA, CB, CG, CG_CD2_length, CB_CG_CD2_angle, CA_CB_CG_CD2_diangle);
	CD2 = Atom(m_atomid+7, "CD2", m_resname, m_resid, cd2);
	Vector3d ce1 = calCoordinates(CB, CG, CD1, CD1_CE1_length, CG_CD1_CE1_angle, CB_CG_CD1_CE1_diangle);
	CE1 = Atom(m_atomid+8, "CE1", m_resname, m_resid, ce1);
	Vector3d ce2 = calCoordinates(CB, CG, CD2, CD2_CE2_length, CG_CD2_CE2_angle, CB_CG_CD2_CE2_diangle);
	CE2 = Atom(m_atomid+9, "CE2", m_resname, m_resid, ce2);
	Vector3d cz = calCoordinates(CG, CD1, CE1, CE1_CZ_length, CD1_CE1_CZ_angle, CG_CD1_CE1_CZ_diangle);
	CZ = Atom(m_atomid+10, "CZ", m_resname, m_resid, cz);
	Vector3d oh = calCoordinates(CD1, CE1, CZ, CZ_OH_length, CE1_CZ_OH_angle, CD1_CE1_CZ_OH_diangle);
	OH = Atom(m_atomid+11, "OH", m_resname, m_resid, oh);
}
void TyrGeo::storeSideChain(RotamerLib& rotamerlib){
	rotamerlib.m_atoms.push_back(CG);
	rotamerlib.m_atoms.push_back(CD1);
	rotamerlib.m_atoms.push_back(CD2);
	rotamerlib.m_atoms.push_back(CE1);
	rotamerlib.m_atoms.push_back(CE2);
	rotamerlib.m_atoms.push_back(CZ);
	rotamerlib.m_atoms.push_back(OH);
}
void TyrGeo::output(vector<Atom>& atomsData){
	GeoContainCB::output(atomsData);
	atomsData.push_back(CG);
	atomsData.push_back(CD1);
	atomsData.push_back(CD2);
	atomsData.push_back(CE1);
	atomsData.push_back(CE2);
	atomsData.push_back(CZ);
	atomsData.push_back(OH);
}

TrpGeo::TrpGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 110.8914;

	C_O_length = 1.23;
	CA_C_O_angle = 120.5117;
	N_CA_C_O_diangle = 120.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	peptide_bond = 1.33;
	CA_C_N_angle = 116.642992978143;
	C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.6112;

	CB_CG_length = 1.50;
	CA_CB_CG_angle = 114.10;
	N_CA_CB_CG_diangle = -66.4;

	CG_CD1_length = 1.37;
	CB_CG_CD1_angle = 127.07;
	CA_CB_CG_CD1_diangle = 96.3;

	CG_CD2_length = 1.43;
	CB_CG_CD2_angle = 126.66;
	CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle-180.0;

	CD1_NE1_length = 1.38;
	CG_CD1_NE1_angle = 108.5;
	CB_CG_CD1_NE1_diangle = 180.0;

	CD2_CE2_length = 1.40;
	CG_CD2_CE2_angle = 108.5;
	CB_CG_CD2_CE2_diangle = 180.0;

	CD2_CE3_length = 1.40;
	CG_CD2_CE3_angle = 133.83;
	CB_CG_CD2_CE3_diangle = 0.0;

	CE2_CZ2_length = 1.40;
	CD2_CE2_CZ2_angle = 120.0;
	CG_CD2_CE2_CZ2_diangle = 180.0;

	CE3_CZ3_length = 1.40;
	CD2_CE3_CZ3_angle = 120.0;
	CG_CD2_CE3_CZ3_diangle = 180.0;

	CZ2_CH2_length = 1.40;
	CE2_CZ2_CH2_angle = 120.0;
	CD2_CE2_CZ2_CH2_diangle = 0.0;

	m_resname =  'W';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 5.4;
	m_num_dasfs = 2;
}
void TrpGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG_diangle = rotamers[0];
		CA_CB_CG_CD1_diangle = rotamers[1];
		if (CA_CB_CG_CD1_diangle > 0){
			CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle-180.0;
		}else{
			CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle+180.0;
		}
	}catch(...){
		cout << "TrpGeo rotamers list: not long enough" << endl;
	}
}
void TrpGeo::addSideChain(){
	Vector3d cg = calCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
	CG = Atom(m_atomid+5, "CG", m_resname, m_resid, cg);
	Vector3d cd1 = calCoordinates(CA, CB, CG, CG_CD1_length, CB_CG_CD1_angle, CA_CB_CG_CD1_diangle);
	CD1 = Atom(m_atomid+6, "CD1", m_resname, m_resid, cd1);
	Vector3d cd2 = calCoordinates(CA, CB, CG, CG_CD2_length, CB_CG_CD2_angle, CA_CB_CG_CD2_diangle);
	CD2 = Atom(m_atomid+7, "CD2", m_resname, m_resid, cd2);
	Vector3d ne1 = calCoordinates(CB, CG, CD1, CD1_NE1_length, CG_CD1_NE1_angle, CB_CG_CD1_NE1_diangle);
	NE1 = Atom(m_atomid+8, "NE1", m_resname, m_resid, ne1);
	Vector3d ce2 = calCoordinates(CB, CG, CD2, CD2_CE2_length, CG_CD2_CE2_angle, CB_CG_CD2_CE2_diangle);
	CE2 = Atom(m_atomid+9, "CE2", m_resname, m_resid, ce2);
	Vector3d ce3 = calCoordinates(CB, CG, CD2, CD2_CE3_length, CG_CD2_CE3_angle, CB_CG_CD2_CE3_diangle);
	CE3 = Atom(m_atomid+10, "CE3", m_resname, m_resid, ce3);
	Vector3d cz2 = calCoordinates(CG, CD2, CE2, CE2_CZ2_length, CD2_CE2_CZ2_angle, CG_CD2_CE2_CZ2_diangle);
	CZ2 = Atom(m_atomid+11, "CZ2", m_resname, m_resid, cz2);
	Vector3d cz3 = calCoordinates(CG, CD2, CE3, CE3_CZ3_length, CD2_CE3_CZ3_angle, CG_CD2_CE3_CZ3_diangle);
	CZ3 = Atom(m_atomid+12, "CZ3", m_resname, m_resid, cz3);
	Vector3d ch2 = calCoordinates(CD2, CE2, CZ2, CZ2_CH2_length, CE2_CZ2_CH2_angle, CD2_CE2_CZ2_CH2_diangle);
	CH2 = Atom(m_atomid+13, "CH2", m_resname, m_resid, ch2);
}
void TrpGeo::storeSideChain(RotamerLib& rotamerlib){
	rotamerlib.m_atoms.push_back(CG);
	rotamerlib.m_atoms.push_back(CD1);
	rotamerlib.m_atoms.push_back(CD2);
	rotamerlib.m_atoms.push_back(NE1);
	rotamerlib.m_atoms.push_back(CE2);
	rotamerlib.m_atoms.push_back(CE3);
	rotamerlib.m_atoms.push_back(CZ2);
	rotamerlib.m_atoms.push_back(CZ3);
	rotamerlib.m_atoms.push_back(CH2);
}
void TrpGeo::output(vector<Atom>& atomsData){
	GeoContainCB::output(atomsData);
	atomsData.push_back(CG);
	atomsData.push_back(CD1);
	atomsData.push_back(CD2);
	atomsData.push_back(NE1);
	atomsData.push_back(CE2);
	atomsData.push_back(CE3);
	atomsData.push_back(CZ2);
	atomsData.push_back(CZ3);
	atomsData.push_back(CH2);
}

