#include <iostream> 
#include <Eigen/Dense>
#include <string>
#include "atom.h"

using namespace std;
using namespace Eigen;

Atom::Atom(){}

Atom::Atom(int param_id, const string& param_name, char param_resname, int param_resid, const Vector3d& param_position){
	m_id = param_id;
	m_name = param_name;
	m_resname = param_resname;
	m_resid = param_resid;
	m_position = param_position;
	if(m_name.compare("N") == 0 || m_name.compare("CA") == 0 || m_name.compare("C") == 0 || m_name.compare("O") == 0){
		m_isMainChain = true;
	}else{
		m_isMainChain = false;
	}
	m_LJType = -1;
	m_radii = -1;
	m_well = -1;

	setLJParams();
}

void Atom::setLJParams(){
	if(m_name.compare("N") == 0){
        if(m_resname == 'P'){
            this->m_LJType = 12;
		}else{
            this->m_LJType = 9;
		}
    }else if(m_name.compare("CA") == 0){
        this->m_LJType = 1;
    }else if(m_name.compare("C") == 0){
        this->m_LJType = 2;
    }else if(m_name.compare("O") == 0){
        this->m_LJType = 13;
    }else if(m_resname == 'A'){
        if (m_name.compare("CB") == 0){
            this->m_LJType = 5;
		}
    }else if(m_resname == 'V'){
        if (m_name.compare("CB") == 0){
            this->m_LJType = 3;
        }else if(m_name.compare("CG1") == 0){
            this->m_LJType = 5;
        }else if(m_name.compare("CG2") == 0){
            this->m_LJType = 5;
		}
    }else if(m_resname == 'I'){
        if (m_name.compare("CB") == 0){
            this->m_LJType = 3;
        }else if(m_name.compare("CG1") == 0){
            this->m_LJType = 5;
        }else if(m_name.compare("CG2") == 0){
            this->m_LJType = 4;
        }else if(m_name.compare("CD1") == 0 || m_name.compare("CD") == 0){
            this->m_LJType = 5;
		}
    }else if(m_resname == 'L'){
        if (m_name.compare("CB") == 0){
            this->m_LJType = 4;
        }else if(m_name.compare("CG") == 0){
            this->m_LJType = 3;
        }else if(m_name.compare("CD1") == 0){
            this->m_LJType = 5;
        }else if(m_name.compare("CD2") == 0){
            this->m_LJType = 5;
		}
    }else if(m_resname == 'S'){
        if (m_name.compare("CB") == 0){
            this->m_LJType = 4;
        }else if(m_name.compare("OG") == 0){
            this->m_LJType = 16;
		}
    }else if(m_resname == 'T'){
        if (m_name.compare("CB") == 0){
            this->m_LJType = 3;
        }else if(m_name.compare("OG1") == 0){
            this->m_LJType = 16;
        }else if(m_name.compare("CG2") == 0){
            this->m_LJType = 5;
		}
    }else if(m_resname == 'D'){
        if (m_name.compare("CB") == 0){
            this->m_LJType = 4;
        }else if(m_name.compare("CG") == 0){
            this->m_LJType = 7;
        }else if(m_name.compare("OD1") == 0){
            this->m_LJType = 15;
        }else if(m_name.compare("OD2") == 0){
            this->m_LJType = 16;
		}
    }else if(m_resname == 'N'){
        if (m_name.compare("CB") == 0){
            this->m_LJType = 4;
        }else if(m_name.compare("CG") == 0){
            this->m_LJType = 7;
        }else if(m_name.compare("OD1") == 0){
            this->m_LJType = 14;
        }else if(m_name.compare("ND2") == 0){
            this->m_LJType = 10;
		}
    }else if(m_resname == 'E'){
        if (m_name.compare("CB") == 0){
            this->m_LJType = 4;
        }else if(m_name.compare("CG") == 0){
            this->m_LJType = 4;
        }else if(m_name.compare("CD") == 0){
            this->m_LJType = 7;
        }else if(m_name.compare("OE1") == 0){
            this->m_LJType = 15;
        }else if(m_name.compare("OE2") == 0){
            this->m_LJType = 16;
		}
    }else if(m_resname == 'Q'){
        if (m_name.compare("CB") == 0){
            this->m_LJType = 4;
        }else if(m_name.compare("CG") == 0){
            this->m_LJType = 4;
        }else if(m_name.compare("CD") == 0){
            this->m_LJType = 7;
        }else if(m_name.compare("OE1") == 0){
            this->m_LJType = 14;
        }else if(m_name.compare("NE2") == 0){
            this->m_LJType = 10;
		}
    }else if(m_resname == 'K'){
        if (m_name.compare("CB") == 0){
            this->m_LJType = 4;
        }else if(m_name.compare("CG") == 0){
            this->m_LJType = 4;
        }else if(m_name.compare("CD") == 0){
            this->m_LJType = 4;
        }else if(m_name.compare("CE") == 0){
            this->m_LJType = 4;
        }else if(m_name.compare("NZ") == 0){
            this->m_LJType = 10;
		}
    }else if(m_resname == 'R'){
        if (m_name.compare("CB") == 0){
            this->m_LJType = 4;
        }else if(m_name.compare("CG") == 0){
            this->m_LJType = 4;
        }else if(m_name.compare("CD") == 0){
            this->m_LJType = 4;
        }else if(m_name.compare("NE") == 0){
            this->m_LJType = 10;
        }else if(m_name.compare("CZ") == 0){
            this->m_LJType = 7;
        }else if(m_name.compare("NH1") == 0){
            this->m_LJType = 10;
        }else if(m_name.compare("NH2") == 0){
            this->m_LJType = 10;
		}
    }else if(m_resname == 'C'){
        if (m_name.compare("CB") == 0){
            this->m_LJType = 8;
        }else if(m_name.compare("SG") == 0){
            this->m_LJType = 17;
		}
    }else if(m_resname == 'M'){
        if (m_name.compare("CB") == 0){
            this->m_LJType = 4;
        }else if(m_name.compare("CG") == 0){
            this->m_LJType = 4;
        }else if(m_name.compare("SD") == 0){
            this->m_LJType = 18;
        }else if(m_name.compare("CE") == 0){
            this->m_LJType = 5;
		}
    }else if(m_resname == 'F'){
        if (m_name.compare("CB") == 0){
            this->m_LJType = 4;
        }else if(m_name.compare("CG") == 0){
            this->m_LJType = 6;
        }else if(m_name.compare("CD1") == 0){
            this->m_LJType = 6;
        }else if(m_name.compare("CD2") == 0){
            this->m_LJType = 6;
        }else if(m_name.compare("CE1") == 0){
            this->m_LJType = 6;
        }else if(m_name.compare("CE2") == 0){
            this->m_LJType = 6;
        }else if(m_name.compare("CZ") == 0){
            this->m_LJType = 6;
		}
    }else if(m_resname == 'Y'){
        if (m_name.compare("CB") == 0){
            this->m_LJType = 4;
        }else if(m_name.compare("CG") == 0){
            this->m_LJType = 6;
        }else if(m_name.compare("CD1") == 0){
            this->m_LJType = 6;
        }else if(m_name.compare("CD2") == 0){
            this->m_LJType = 6;
        }else if(m_name.compare("CE1") == 0){
            this->m_LJType = 6;
        }else if(m_name.compare("CE2") == 0){
            this->m_LJType = 6;
		}else if (m_name.compare("CZ") == 0){
            this->m_LJType = 6;
		}else if (m_name.compare("OH") == 0){
            this->m_LJType = 16;
		}
    }else if(m_resname == 'W'){
        if (m_name.compare("CB") == 0){
            this->m_LJType = 4;
        }else if(m_name.compare("CG") == 0){
            this->m_LJType = 6;
        }else if(m_name.compare("CD1") == 0){
            this->m_LJType = 6;
        }else if(m_name.compare("CD2") == 0){
            this->m_LJType = 6;
        }else if(m_name.compare("NE1") == 0){
            this->m_LJType = 11;
        }else if(m_name.compare("CE2") == 0){
            this->m_LJType = 6;
        }else if(m_name.compare("CE3") == 0){
            this->m_LJType = 6;
        }else if(m_name.compare("CZ2") == 0){
            this->m_LJType = 6;
        }else if(m_name.compare("CZ3") == 0){
            this->m_LJType = 6;
        }else if(m_name.compare("CH2") == 0){
            this->m_LJType = 6;
		}
    }else if(m_resname == 'H'){
        if (m_name.compare("CB") == 0){
            this->m_LJType = 4;
        }else if(m_name.compare("CG") == 0){
            this->m_LJType = 6;
        }else if(m_name.compare("ND1") == 0){
            this->m_LJType = 11;
        }else if(m_name.compare("CD2") == 0){
            this->m_LJType = 6;
        }else if(m_name.compare("CE1") == 0){
            this->m_LJType = 6;
        }else if(m_name.compare("NE2") == 0){
            this->m_LJType = 11;
		}
    }else if(m_resname == 'P'){
        if (m_name.compare("CB") == 0){
            this->m_LJType = 6;
        }else if(m_name.compare("CG") == 0){
            this->m_LJType = 6;
        }else if(m_name.compare("CD") == 0){
            this->m_LJType = 6;
		}
	}

	if (this->m_LJType == -1 && m_name[0] != 'H'){
		cout << "LJ params not found: " + to_string(m_resid) + " " + m_resname + " " + m_name << endl;
	}else{
		this->m_radii = RADII[this->m_LJType-1];
		this->m_well = WELL[this->m_LJType-1];
	}

}





/*
int main(){
	Vector3d v(1,1,1);
	Atom a(10,"CB",'A',10,v);

	cout << a.m_id << endl;
	cout << a.m_name << endl;
	cout << a.m_resname << endl;
	cout << a.m_LJType << endl;
	cout << a.m_radii << endl;
	cout << a.m_well << endl;
	cout << a.m_isMainChain << endl;

	a.setLJParams('M',"SD");
	
	cout << RADII[17] << endl;
	cout << a.m_LJType << endl;
	cout << a.m_radii << endl;
	cout << a.m_well << endl;
	return 0;
}*/


