#ifndef PEPTIDEBUILDER_H
#define PEPTIDEBUILDER_H

#include <iostream>
#include <vector> 
#include "residue.h"
#include "atom.h"
#include "geometry.h"

using namespace std;

vector<Residue> initResidueChain(const vector<char>& resnames);

Geo* getGeo(int resid, char resname, int* totalAtoms);

vector<Geo*> initGeoChain(vector<Residue>& residuesData);

vector<Atom> rebuildSideChain(vector<Geo*>& geosData, vector<Residue>& residuesData);

vector<Atom> rebuildMainChainCb(vector<Geo*>& geosData, vector<Residue>& residuesData);

vector<Atom> outputAtomsData(const vector<Geo*>& geosData);

#endif