#ifndef PDBIO_H
#define PDBIO_H

#include<string> 
#include <vector> 
#include "atom.h"

using namespace std;

void outputPDB(const string& outfile, const vector<Atom>& atomsData);

vector<Atom> readPDB(const string& infile);

#endif

