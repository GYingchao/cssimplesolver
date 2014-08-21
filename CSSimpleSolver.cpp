#include "CSSimpleSolver.h"
using namespace std;

bool CSSimpleSolver::parser(FILE *f, vector<vector<double>> &A, vector<double> &b)
{
	if(f == NULL) return false;
	// Read the input data and allocate memory
	// ...
	return true;
}

bool CSSimpleSolver::L1Solver(vector<vector<double>> &A, vector<double> &b)
{
	// Find a proper left null space of A : F
	//-- [Assume that F is sparse...] --//
	vector<vector<double>> F = leftNullSpace(A);
	return false;
}

bool CSSimpleSolver::OMPSolver(vector<vector<double>> &A, vector<double> &b)
{
	return false;
}

// Compute Ax = 0, return basis vectors of solution space 
// A is assumed to be sparse matrix
vector<vector<double>>& CSSimpleSolver::nullSpace(const vector<vector<double>> &A)
{
	vector<vector<double>> tem;
	if(!A.empty()) {
		int m = A.size();
		int n = A[0].size();
		/*
		//if(m >= n) {
			// Can be solved by cholmod
			Cholmod_Solver solver;
			SparseMatrix S(m, n, false);
			DenseMatrix x(n, 1);
			DenseMatrix o(m, 1);
			for(int i=0; i<m; i++) {
				o.set(i, 0, 0.0);
				for(int j=0; j<n; j++) {
					S.add_coef(i, j, A[i][j]);
				}
			}
			// Solve Sx = o
			if(!solver.linear_solver(S, o, x)) {
				cout << "Error: Solving nullspace failed!" << endl;
			} else {
				// Naively, construct the nullspace by basic unit vectors of x
				int xm = x.size();
				int xn = x[0].size();
				for(int i=0; i<xm; i++) {
					for(int j=0; j<xn; j++) {
						tem[i][j] = x.get(i, j);
					}
				}
			}
		//} else {
			//cout << "Matrix format error!" << endl;
		//}
		*/
	} else {
		cout << "Null pointer of A in computing nullspace!" << endl;
	}
	return tem;
}

// Compute FA = 0, return proper F.
vector<vector<double>>& CSSimpleSolver::leftNullSpace(const vector<vector<double>> &A)
{
	/*
		Naively:
		To compute F s.t. FA = 0, by Associative Law of Matrix Transpose, A'F' = 0;
		We just need to compute the nullspace of A' to get F', then F.
	*/
	vector<vector<double>> F;
	vector<vector<double>> A_T;
	int m = A.size();
	int n = A[0].size();
	for(int i=0; i<n; i++) {
		for(int j=0; j<m; j++) {
			A_T[i][j] = A[j][i];
		}
	}

	vector<vector<double>> x = nullSpace(A_T);

	if(!x.empty()) {
		int xm = x.size();
		int xn = x[0].size();
		for(int i=0; i<xn; i++) {
			for(int j=0; j<xm; j++) {
				F[i][j] = x[j][i];
			}
		}
	} else {
		cout << "Error: Compute left null space failed!" << endl;
	}
	return F;
}


//--- For test ---//
int main()
{
	// Test nullSpace()
	vector<vector<double>> A;
	double a0[] = {-1, 1, 2, 4};
	double a1[] = {2, 0, 1, -7};
	vector<double> A0(a0, a0 + sizeof(a0)/sizeof(double));
	vector<double> A1(a1, a1 + sizeof(a1)/sizeof(double));
	A.push_back(A0); A.push_back(A1);

	CSSimpleSolver css;
	vector<vector<double>> x = css.nullSpace(A);
	for(int i=0; i<x.size(); i++) {
		for(int j=0; j<x[0].size()-1; j++) {
			cout << x[i][j] << "\t" ;
		}
		cout << x[i][x[0].size()-1] << endl;
	}
	return 0;
}
