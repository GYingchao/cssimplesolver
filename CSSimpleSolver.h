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
#include "lp_lib.h"

#define min(a,b) ((a)>(b)?(b):(a))
#define max(a,b) ((a)<(b)?(b):(a))

using namespace std;

#define CPUTIME (clock())

/*
typedef Cholmod_solver_traits<double> Cholmod_Solver;
typedef Cholmod_Solver::Dense_matrix DenseMatrix;
typedef Cholmod_Solver::Sparse_matrix SparseMatrix;
*/


/*
Solver for the Following Compressed Sensing Model:
Ax + e = b, where A(m, n) : m>n; e is sparse.
*/
class CSSimpleSolver
{
public:
	CSSimpleSolver() {}
	~CSSimpleSolver() {}

	bool parser(FILE *f, vector<vector<double>> &A, vector<double> &b);

	bool L1Solver(vector<vector<double>> &A, vector<double> &b);
	bool OMPSolver(vector<vector<double>> &A, vector<double> &b);

	bool LPSolver(vector<vector<double>> &A, vector<double> &b, vector<double> &CT, vector<double> &x, double &optimal);

	bool Solve(vector<vector<double>> &A, vector<double> &b, vector<double> &x, vector<double>& e);

//private:
	//cholmod_common c;
	vector<vector<double>> nullSpace(const vector<vector<double>> &A);
	vector<vector<double>> leftNullSpace(const vector<vector<double>> &A);
	void print_matrix( char* desc, int m, int n, double* a, int lda );
	bool multiply_matrix(vector<vector<double>> &A, vector<vector<double>> &B, vector<vector<double>> &C);
};
#endif