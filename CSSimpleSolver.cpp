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
vector<vector<double>> CSSimpleSolver::nullSpace(const vector<vector<double>> &A)
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

		/* Use lapacke to compute nullspace by calling SVD */
		double* a = new double[m*n];
		int aIndex = 0;
		for(int i=0; i<m; i++) {
			for(int j=0; j<n; j++) {
				a[aIndex++] = A[i][j];
			}
		}

		double* s = new double[n];
		double* u = new double[m*m];
		double* vt = new double[n*n];
		double* superb = new double[min(m, n) -1];
		int lda = max(m ,n), ldu = m, ldvt = n, info;

		// Compute SVD
		info = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, 'A', 'A', m, n, a, lda,
                        s, u, ldu, vt, ldvt, superb);

		// Check result
		if(info > 0) {
			printf("Error in computing nullspace: SVD failed to converge!\n");
			exit(-1);
		}

		print_matrix("U: ", m, n, u, m);
		print_matrix("Singular Values: ", 1, n, s, 1);
		print_matrix("VT: ", n, n, vt, n);

		// Define a very small value eps, and set the lower bound to eps*min(m, n), which is similar to Matlab implementation
		//double eps = std::numeric_limits<float>::denorm_min()*max(l_m, l_n);
		double eps = 0.00000000001;
		printf("eps is %1.11f\n", eps);
		// Check the small enough singular values to get the corresponding nullspace vectors
		int startIndex = 0;
		for(int i=n-1; i>=0; i--) {
			if(s[i] - eps < eps) continue;
			else {
				startIndex = i+1;
				break;
			}
		}

		/*for (int i = 0; i < l_n; i++) {
			if (s[i] > eps) continue;
			else {
				startIndex = i;
			}
		}*/

		/*int vStart = startIndex*n;
		int vEnd = n*n;
		double* null;
		for(int i = vStart; i<vEnd; i++) {
			null[i-vStart] = vt[i];
		}*/

		delete s, u, vt;
		printf("startIndex %d\n", startIndex);

		// Construct nullspace matrix
		int vIndex = startIndex*n;
		//for(int i=0; i<n-startIndex; i++) {
		for(int i=0; i<n; i++) {
			tem.push_back(vector<double>());
			//for(int j=0; j<n; j++) {
			for(int j=0; j<n-startIndex; j++) {
				tem[i].push_back(vt[vIndex++]);
			}
		}
		
	} else {
		cout << "Null pointer of A in computing nullspace!" << endl;
	}
	cout << "tem size: " << tem.size() << endl;
	return tem;
}

// Compute FA = 0, return proper F.
vector<vector<double>> CSSimpleSolver::leftNullSpace(const vector<vector<double>> &A)
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
		A_T.push_back(vector<double>());
		for(int j=0; j<m; j++) {
			A_T[i].push_back(A[j][i]);
		}
	}

	vector<vector<double>> x = nullSpace(A_T);

	if(!x.empty()) {
		int xm = x.size();
		int xn = x[0].size();
		for(int i=0; i<xn; i++) {
			F.push_back(vector<double>());
			for(int j=0; j<xm; j++) {
				F[i].push_back(x[j][i]);
			}
		}
	} else {
		cout << "Error: Compute left null space failed!" << endl;
	}
	return F;
}

void CSSimpleSolver::print_matrix( char* desc, int m, int n, double* a, int lda ) {
        lapack_int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.8f", a[i*lda+j] );
                printf( "\n" );
        }
}
//--- For test ---//
int main()
{
	//// Test nullSpace()
	//vector<vector<double>> A;
	//double a0[] = {-1, 1, 2, 4};
	//double a1[] = {2, 0, 1, -7};
	//vector<double> A0(a0, a0 + sizeof(a0)/sizeof(double));
	//vector<double> A1(a1, a1 + sizeof(a1)/sizeof(double));
	//A.push_back(A0); A.push_back(A1);

	// Test leftNullSpace()
	vector<vector<double>> A;
	double a0[] = {-1, 2};
	double a1[] = {1, 0};
	double a2[] = {2, 1};
	double a3[] = {4, -7};
	vector<double> A0(a0, a0 + sizeof(a0)/sizeof(double));
	vector<double> A1(a1, a1 + sizeof(a1)/sizeof(double));
	vector<double> A2(a2, a2 + sizeof(a2)/sizeof(double));
	vector<double> A3(a3, a3 + sizeof(a3)/sizeof(double));
	A.push_back(A0); A.push_back(A1); A.push_back(A2); A.push_back(A3);


	CSSimpleSolver css;
	vector<vector<double>> x = css.leftNullSpace(A);
	for(int i=0; i<x.size(); i++) {
		for(int j=0; j<x[0].size()-1; j++) {
			cout << x[i][j] << "\t" ;
		}
		cout << x[i][x[0].size()-1] << endl;
	}
	return 0;
}
