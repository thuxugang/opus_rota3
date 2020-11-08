#ifndef ATOM_H
#define ATOM_H

#include<string> 
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

const double RADII[18] = {2, 2, 2, 2, 2, 1.8, 2, 2, 1.75, 1.75, 1.75, 1.75, 1.55, 1.55, 1.55, 1.55, 1.9, 1.9};
const double WELL[18] = {0.0486, 0.14, 0.0486, 0.1142, 0.1811, 0.08, 0.12, 0.1142, 0.2384, 0.2384, 0.2384, 0.2384, 0.1591, 0.1591, 0.21, 0.1591, 0.16, 0.16};

class Atom{
public:
	int m_id;
	string m_name;
	char m_resname;
	int m_resid;
	char m_restype;
	Vector3d m_position;

	int m_LJType;
	double m_radii;
	double m_well;

	bool m_isMainChain;

	Atom();
	Atom(int param_id, const string& param_name, char param_resname, int param_resid, const Vector3d& param_position);
	void setLJParams();
};

#endif