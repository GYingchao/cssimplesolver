#ifndef CSSIMPLESOLVER_H_
#define CSSIMPLESOLVER_H_

#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
//#include <Cholmod_solver_traits.h>
#include "lapacke.h"

#define min(a,b) ((a)>(b)?(b):(a))
#define max(a,b) ((a)<(b)?(b):(a))

using namespace std;

#define CPUTIME (clock())

/*
typedef Cholmod_solver_traits<double> Cholmod_Solver;
typedef Cholmod_Solver::Dense_matrix DenseMatrix;
typedef Cholmod_Solver::Sparse_matrix SparseMatrix;
*/

class CSSimpleSolver
{
public:
	CSSimpleSolver() {}
	~CSSimpleSolver() {}

	bool parser(FILE *f, vector<vector<double>> &A, vector<double> &b);

	bool L1Solver(vector<vector<double>> &A, vector<double> &b);
	bool OMPSolver(vector<vector<double>> &A, vector<double> &b);

//private:
	//cholmod_common c;
	vector<vector<double>>& nullSpace(const vector<vector<double>> &A);
	vector<vector<double>>& leftNullSpace(const vector<vector<double>> &A);
};
#endif