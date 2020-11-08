#include "residue.h"
#include "rotamerLibScore.h"

using namespace std;

double getRLScore(Residue& residue, int index){
	if(residue.m_resname == 'G' || residue.m_resname == 'A'){
		return 0.0;
	}else{
		double p = residue.m_rotamerlibs[index].m_prob;
		double p0 = residue.m_rotamerlibs[0].m_prob;
		double rlscore = -log(p/p0);
		rlscore = rlscore>5?5:rlscore;
		rlscore = rlscore<-5?-5:rlscore;
		return rlscore;
	}
}