#ifndef CONFIG_H
#define CONFIG_H

#include <map>  
#include <string> 
#include <random>

using namespace std;

extern map<string, string> CONFIG;
extern default_random_engine E;

vector<string> split(string& s, const string& c);
string trim(const string& str);
void readConfigFile(string path);

#endif

