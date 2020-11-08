#include "otherRL.h"

vector<string> splits(stringstream& ss, const string &s, char delim){
	ss << s;
	string item;
	vector<string> elems;
	while (getline(ss, item, delim)){
		elems.push_back(item);
	}
	ss.clear();
	return elems;
}


void addRotamersInfoFromFASPR(const string& inputname, const string& tmp_dir, vector<Residue>& residuesData){

	string outputfile = tmp_dir + "/faspr.chi";
	system(("./FASPR/FASPR -i " + inputname + " -o " + outputfile).c_str());

	cout << "FASPR Done..." << endl;

	stringstream ss;

	ifstream fasprFile;
	try{
		fasprFile.open(outputfile);
	}
	catch (...){
		cout << "FASPR output error!";
		fasprFile.close();
		exit(-1);
	}

	map<int, Residue*> m_residue_dict;
	for (Residue& r : residuesData){
		m_residue_dict[r.m_resid] = &r;
	}

	string str_line;
	try{
		while (!fasprFile.eof()){

			getline(fasprFile, str_line);

			if (str_line.find('#') == 0 || str_line.find('=') == 0 || str_line.size() == 0){
				continue;
			}

			str_line = trim(str_line);

			char spliter = ' ';
			if (str_line.find('\t') != str_line.npos){
				spliter = '\t';
			}

			vector<string> items = splits(ss, str_line, spliter);

			Residue* r = m_residue_dict[atof(items[0].c_str())];

			//assert(r->m_resname != getResname(items[1]));
			if (r->m_resname != getResname(items[1])){
				throw 0;
			}

			vector<float> chi = { 0, 0, 0, 0 };
			for (int i = 2; i < items.size(); ++i){
				chi[i - 2] = atof(items[i].c_str());
			}
			RotamerLib rl(-2, stof(CONFIG["w_faspr"])*r->m_rotamerlibs[r->fisrt_rl_id].m_prob, chi[0], chi[1], chi[2], chi[3]);

			//insert at first
			r->m_rotamerlibs.insert(r->m_rotamerlibs.begin(), rl);
			r->fisrt_rl_id += 1;

		}
	}catch (...){
		cout << "FASPR format error!" << endl;
		fasprFile.close();
	}
	fasprFile.close();

}

void calOSCAR(const string& list_path){
	system(("cp " + list_path + " ./oscar/data && \
		cd ./oscar && ./oscar-star > oscar.log").c_str());
}

void addRotamersInfoFromOSCAR(string inputname, vector<Residue>& residuesData){

	vector<string> names = split(inputname, "/");
	string name = names[names.size() - 1];
	
	names = split(name, ".");
	name = names[0];

	//system(("cd ./oscar && \
	//	./oscar-star " + inputname + " > oscar.log").c_str());

	//cout << "OSCAR-star Done..." << endl;

	stringstream ss;

	vector<Atom> ad;
	ifstream oscarFile;
	try{
		ad = readPDB("./oscar/" + name + "_model.pdb");
	}catch (...){
		cout << "OSCAR output error!";
		oscarFile.close();
		exit(-1);
	}

	system(("rm ./oscar/" + name + "_model.pdb").c_str());

	vector<Residue> rd = atomToResidue(ad, false);
	addDihedral(rd);
	assert(rd.size() == residuesData.size());
	
	try{
		for (int i = 0; i < rd.size(); ++i){

			if (residuesData[i].m_resname == 'G' || residuesData[i].m_resname == 'A') continue;

			vector<float> chi = { 0, 0, 0, 0 };
			for (int j = 0; j < rd[i].m_dihedrals.size(); ++j){
				chi[j] = rd[i].m_dihedrals[j];
			}

			RotamerLib rl(-4, stof(CONFIG["w_oscar"])*residuesData[i].m_rotamerlibs[residuesData[i].fisrt_rl_id].m_prob, chi[0], chi[1], chi[2], chi[3]);

			//insert at first
			residuesData[i].m_rotamerlibs.insert(residuesData[i].m_rotamerlibs.begin(), rl);
			residuesData[i].fisrt_rl_id += 1;

		}
	}catch (...){
		cout << "OSCAR format error!" << endl;
	}
}

map<string, string> getConstraintList(string& path){

	map<string, string> list;

	ifstream file;
	file.open(path.c_str());
	string str_line;
	try{
		while (!file.eof()){
			getline(file, str_line);
			str_line = trim(str_line);
			size_t pos = str_line.find(' ');
			string str_key = str_line.substr(0, pos);
			string str_value = str_line.substr(pos + 1);
			list.insert(pair<string, string>(trim(str_key), trim(str_value)));
		}
	}catch (...){
		cout << path << ": constraint list format error!" << endl;
		file.close();
		exit(-1);
	}
	file.close();

	return list;
}

void addRotamersInfoFromOthers(const string& inputname, const map<string, string>& other_list, vector<Residue>& residuesData){

	stringstream ss;

	ifstream othersFile;
	try{
		othersFile.open(other_list.at(inputname));
	}catch (...){
		cout << inputname << ": format error!";
		othersFile.close();
		exit(-1);
	}

	map<int, Residue*> m_residue_dict;
	for (Residue& r : residuesData){
		m_residue_dict[r.m_resid] = &r;
	}

	string str_line;
	try{
		while (!othersFile.eof()){

			getline(othersFile, str_line);

			if (str_line.find('#') == 0 || str_line.find('=') == 0 || str_line.size() == 0){
				continue;
			}

			str_line = trim(str_line);

			char spliter = ' ';
			if (str_line.find('\t') != str_line.npos){
				spliter = '\t';
			}

			vector<string> items = splits(ss, str_line, spliter);
			assert(items.size() == 6);

			Residue* r = m_residue_dict[atof(items[0].c_str())];

			if (r->m_resname == 'G' || r->m_resname == 'A') continue;

			//assert(r->m_resname != getResname(items[1]));
			if (r->m_resname != getResname(items[1])){
				throw 0;
			}

			vector<float> chi = { 0, 0, 0, 0 };
			for (int i = 2; i < items.size(); ++i){
				if (atof(items[i].c_str()) <= 180){
					chi[i - 2] = atof(items[i].c_str());
				}			
			}
			RotamerLib rl(-3, stof(CONFIG["w_others"])*r->m_rotamerlibs[r->fisrt_rl_id].m_prob, chi[0], chi[1], chi[2], chi[3]);

			//insert at first
			r->m_rotamerlibs.insert(r->m_rotamerlibs.begin(), rl);
			r->fisrt_rl_id += 1;
		}
	}catch (...){
		cout << "Others format error!" << endl;
		othersFile.close();
	}
	othersFile.close();

}

void addRotamersInfoFromRotaNN(const string& inputname, const map<string, string>& rotann_list, vector<Residue>& residuesData){

	stringstream ss;

	ifstream rotannFile;
	try{
		rotannFile.open(rotann_list.at(inputname));
	}catch (...){
		cout << inputname << ": format error!";
		rotannFile.close();
		exit(-1);
	}

	vector<string> str_lines;
	string str_line;
	while (!rotannFile.eof()){
		getline(rotannFile, str_line);
		if (str_line.find('#') == 0 || str_line.find('=') == 0 || str_line.size() == 0){
			continue;
		}
		str_lines.push_back(trim(str_line));
	}
	rotannFile.close();

	if (str_lines.size() != residuesData.size()){
		cout << "RotaNN size error!" << endl;
		return;
	}

	try{
		for (int i = 0; i < str_lines.size(); ++i){

			char spliter = ' ';
			if (str_lines[i].find('\t') != str_lines[i].npos){
				spliter = '\t';
			}

			vector<string> items = splits(ss, str_lines[i], spliter);
			assert(items.size() == 5);

			if (residuesData[i].m_resname == 'G' || residuesData[i].m_resname == 'A') continue;

			vector<float> chi = { 0, 0, 0, 0 };
			for (int j = 1; j < items.size(); ++j){
				if (atof(items[j].c_str()) <= 180){
					chi[j - 1] = atof(items[j].c_str());
				}
			}
			RotamerLib rl(-1, stof(CONFIG["w_rotann"])*residuesData[i].m_rotamerlibs[residuesData[i].fisrt_rl_id].m_prob, chi[0], chi[1], chi[2], chi[3]);

			//insert at first
			residuesData[i].m_rotamerlibs.insert(residuesData[i].m_rotamerlibs.begin(), rl);
			residuesData[i].fisrt_rl_id += 1;
		}
	}catch (...){
		cout << "RotaNN format error!" << endl;
	}
}

void addRotamersInfoRefineX4(vector<Residue>& residuesData){
	for (Residue& r : residuesData){

		if (r.m_resname == 'G' || r.m_resname == 'A') continue;

		vector<RotamerLib> tmp_rls;
		for (int i = 0; i <= r.fisrt_rl_id; ++i){
			RotamerLib rl(-5 - i, stof(CONFIG["w_refinex4"])*r.m_rotamerlibs[i].m_prob, r.m_rotamerlibs[i].m_dihedral[0], r.m_rotamerlibs[i].m_dihedral[1],
				r.m_rotamerlibs[i].m_dihedral[2], r.m_rotamerlibs[r.fisrt_rl_id - 1].m_dihedral[3]);
			tmp_rls.push_back(rl);	
		}	
		for (int i = 0; i < tmp_rls.size(); ++i){
			r.m_rotamerlibs.insert(r.m_rotamerlibs.begin() + i, tmp_rls[i]);
		}
	}
}