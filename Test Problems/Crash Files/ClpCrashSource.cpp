#include "Coin_C_defines.h"
#include "ClpSimplex.hpp"

int main()
{
	CoinBigIndex Hcol[3] = {0,1,1};
	int Hrow[1] = {1};
	double H[1] = {1.0};
	double f[2] = {0,0};
	double lb[2] = {-0.5,0.85};
	double ub[2] = {-0.35,1};

	CoinBigIndex Acol[3] = {0,2,4};
	int Arow[4] = {0,1,0,1};
	double A[4] = {1,1,1,-1};
	double rl[2] = {-COIN_DBL_MAX,-COIN_DBL_MAX};
	double ru[2] = {1.5,1.2};

	ClpSimplex simplex;
	simplex.loadProblem(2,2,Acol,Arow,A,lb,ub,f,rl,ru);
	simplex.loadQuadraticObjective(2,Hcol,Hrow,H);
	simplex.setLogLevel(4);
	simplex.initialSolve();
	return 0;
}