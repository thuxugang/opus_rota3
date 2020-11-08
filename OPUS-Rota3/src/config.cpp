#include <fstream> 
#include <iostream> 
#include <map>  
#include <string> 
#include <random>
#include "config.h"

using namespace std;

map<string, string> CONFIG;
default_random_engine E(666);

string lTrim(const string& str) {
	return str.substr(str.find_first_not_of(" \n\r\t"));
}

string rTrim(const string& str) {
	return str.substr(0,str.find_last_not_of(" \n\r\t")+1);
}

string trim(const string& str) {
	return lTrim(rTrim(str));
}

vector<string> split(string& s, const string& c) {

    vector<string> v;
    string::size_type pos1, pos2;
    pos2 = s.find(c);
    pos1 = 0;
    while(string::npos != pos2)
    {
        v.push_back(s.substr(pos1, pos2-pos1));

        pos1 = pos2 + c.size();
        pos2 = s.find(c, pos1);
    }
    if(pos1 != s.length())
        v.push_back(s.substr(pos1));
    return v;
}

void readConfigFile(string path){
	ifstream configFile;
	configFile.open(path.c_str());
	string str_line;
	if (configFile.is_open()){
		while (!configFile.eof()){
			getline(configFile, str_line);
			if (str_line.find('#') == 0 || str_line.find('=') == 0 || str_line == "" || str_line == " " || str_line == "\t" || str_line == "\r"){
				continue;
			}
			size_t pos = str_line.find('=');
			string str_key = str_line.substr(0, pos);
			string str_value = str_line.substr(pos + 1);
			CONFIG.insert(pair<string, string>(trim(str_key), trim(str_value)));
		}
	}else{
		cout << "Cannot open config file" << endl;
		exit(-1);
	}
}