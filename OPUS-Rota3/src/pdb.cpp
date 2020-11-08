#include <fstream> 
#include <iostream> 
#include <sstream>
#include<string> 
#include <vector> 
#include "atom.h"
#include "residue.h"
#include "pdb.h"

using namespace std;

void outputPDB(const string& outfile, const vector<Atom>& atomsData){

	ofstream fout(outfile);

	char data[80];
	int length = atomsData.size();
	for(int i=0; i< length; i++){
		string triresname = getTriResname(atomsData[i].m_resname);
		double x = atomsData[i].m_position[0];
		double y = atomsData[i].m_position[1];
		double z = atomsData[i].m_position[2];
		//sprintf(data, "%-6s%5d%2c%-3s%1c%3s%2c%4d%4c%8.3f%8.3f%8.3f","ATOM",atomsData[i].m_id, ' ', atomsData[i].m_name.c_str(), ' ',triresname.c_str(),' ',atomsData[i].m_resid, ' ',x,y,z);
		sprintf(data, "%-6s%5d%2c%-3s%1c%3s%2c%4d%4c%8.3f%8.3f%8.3f","ATOM",i+1, ' ', atomsData[i].m_name.c_str(), ' ',triresname.c_str(),' ',atomsData[i].m_resid, ' ',x,y,z);
		fout << data << endl;
	}
	fout << flush;
	fout.close();
}

vector<Atom> readPDB(const string& infile){

	vector<Atom> atomsData;

	ifstream fin(infile);
	string current_line;

	stringstream current_val;

	while (!fin.eof()){
				
		getline(fin, current_line);
		//cout << current_line << endl;

		string label = current_line.substr(0, 4);

		if(label.compare("ATOM") == 0){

			int atom_id;
			string atom_name;
			char atom_resname;
			int atom_resid;
			
			char atom_restype;

			current_val.str(current_line.substr(6, 5));
			current_val >> atom_id;
			current_val.clear();

			current_val.str(current_line.substr(12, 4));
			current_val >> atom_name;
			current_val.clear();

			string atom_restype_string = current_line.substr(16, 1);
			if(atom_restype_string.compare(" ") == 0){
				atom_restype = 'A';
			}else{
				current_val.str(atom_restype_string);
				current_val >> atom_restype;
				current_val.clear();
			}

			string triresname;
			current_val.str(current_line.substr(17, 3));
			current_val >> triresname;
			atom_resname = getResname(triresname);
			current_val.clear();

			current_val.str(current_line.substr(22, 4));
			current_val >> atom_resid;
			current_val.clear();

			double x, y ,z;
			current_val.str(current_line.substr(30, 8));
			current_val >> x;
			current_val.clear();
			current_val.str(current_line.substr(38, 8));
			current_val >> y;
			current_val.clear();
			current_val.str(current_line.substr(46, 8));
			current_val >> z;
			current_val.clear();
			Vector3d atom_position(x,y,z);
			
			Atom atom = Atom(atom_id, atom_name, atom_resname, atom_resid, atom_position);
			atom.m_restype = atom_restype;

			atomsData.push_back(atom);
		}		
	}
	fin.close();
	return atomsData;
}

//int main(){
//	//vector<Atom> aa;
//	//Vector3d v(1,1,1);
//	//Atom a(10, "CB", 'G', 1, v);
//	//aa.push_back(a);
//	//outputPDB("test.pdb", aa);
//
//	vector<Atom> ad = readPDB("data/ori1a7s.pdb");
//
//	outputPDB("data/1a7s_rebuild.pdb", ad);
//	return 0;
//
//}
