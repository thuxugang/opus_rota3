#include <iostream> 
#include <string> 
#include <vector> 
#include <map>
#include <time.h>
#include "atom.h"
#include "residue.h"
#include "pdb.h"
#include "angle.h"
#include "peptideBuilder.h"
#include "redis.h"
#include "geometry.h"
#include "potential.h"
#include "config.h"
#include "MonteCarlo.h"
#include "file.h"
#include "otherRL.h"

using namespace std;

int main(){

	cout << "Start..." << endl;

	//read config
	readConfigFile("./rota3.ini");

	cout << "Input dir: " << CONFIG["inputlist"] << endl;
	cout << "Output dir: " << CONFIG["outputdir"] << endl;
	if(stoi(CONFIG["optimized"])){
		cout << "Optimized: true" << endl;
		cout << "Optimization times: " << CONFIG["optimization_times"] << endl;
	}else{
		cout << "Optimized: false" << endl;
	}

	if (stoi(CONFIG["use_rotann"])){
		cout << "Use RotaNN Results..." << endl;
	}
	if (stoi(CONFIG["use_faspr"])){
		cout << "Use FASPR Results..." << endl;
	}
	if (stoi(CONFIG["use_oscar"])){
		cout << "Use OSCAR-star Results..." << endl;
	}
	if (stoi(CONFIG["use_others"])){
		cout << "Use Other Results..." << endl;
	}

	cout << "Connect database..." << endl;
	//connect redis
	Redis r = Redis();
	cout << "Successed!" << endl;

	time_t star_time = time(NULL);

	vector<string> files = getFiles(CONFIG["inputlist"]);
	int length = files.size();
	if (stoi(CONFIG["use_oscar"])){
		cout << "OSCAR-star Start..." << endl;
		calOSCAR(CONFIG["inputlist"]);
		cout << "OSCAR-star Done..." << endl;
	}


	for(int i=0; i< length; i++){
		cout << files[i] << endl;
		string input_file = files[i];
		vector<string> names = split(files[i], "/");
		string output_file = CONFIG["outputdir"]+"/rota3_"+names[names.size()-1];

		//get main chain atom from pdb
		vector<Atom> ad = readPDB(input_file);
		vector<Residue> rd = atomToResidue(ad);
		//init geo
		vector<Geo*> gd = initGeoChain(rd);

		//add phipsi information to residuedata
		addPhiPsi(rd);

		//add rotamer information to residuedata
		addRotamersInfo(rd, r, stoi(CONFIG["top_n"]));

		if (stoi(CONFIG["use_rotann"])){
			map<string, string> rotann_list = getConstraintList(CONFIG["rotann_list"]);
			addRotamersInfoFromRotaNN(input_file, rotann_list, rd);
		}
		if (stoi(CONFIG["use_faspr"])){
			addRotamersInfoFromFASPR(input_file, CONFIG["tmp_dir"], rd);
		}
		if (stoi(CONFIG["use_others"])){
			map<string, string> other_list = getConstraintList(CONFIG["others_list"]);
			addRotamersInfoFromOthers(input_file, other_list, rd);
		}
		if (stoi(CONFIG["use_oscar"])){
			addRotamersInfoFromOSCAR(input_file, rd);
		}

		//add dasf information to residuedata
		float coverage = addDASFInfo(rd, gd, r);

		//get contact list
		vector<Contact> contact_list = getContactList(gd, rd);
		//add all rotamer potentials to residuedata.RotamerLib
		initAllPotentials(gd, rd, coverage);
		//sampling
		getBestStructure(rd);

		//rebuild from m_rl_id
		vector<Atom> ad_rebulid = rebuildSideChain(gd, rd);
		//output
		outputPDB(output_file, ad_rebulid);
		//free Geo*
		int gd_length = gd.size();
		for(int j=0; j<gd_length; j++){
			delete(gd[j]);
		}
		
	}

	time_t end_time = time(NULL);
	cout << "Total time: " << (end_time - star_time) << "s" << endl;

	r.destroy();

	cout << "Done!" << endl;

	return 0;
}