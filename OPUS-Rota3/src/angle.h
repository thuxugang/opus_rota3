#ifndef ANGLE_H
#define ANGLE_H

#include <Eigen/Dense>
#include "residue.h"
#include "atom.h"

using namespace std;
using namespace Eigen;

//const double M_PI = 3.14159265358979323846;

double getNorm(const Vector3d& v);
double getAngle(const Vector3d& v1, const Vector3d& v2);
double calDihedral(const Vector3d& v1, const Vector3d& v2, const Vector3d& v3, const Vector3d& v4);
double calAngle(const Vector3d& v1, const Vector3d& v2, const Vector3d& v3);
MatrixXd rotaxis(double theta, const Vector3d& v);
Vector3d calCoordinates(const Atom& refA, const Atom& refB, const Atom& refC, double L, double ang, double di);

void addPhiPsi(vector<Residue>& residuesData);
void addDihedral(vector<Residue>& residuesData);

#endif