#include "CSSimpleSolver.h"
//using namespace std;

bool CSSimpleSolver::parser(FILE *f, vector<vector<double>> &A, vector<double> &b)
{
	if(f == NULL) return false;
	// Read the input data and allocate memory
	// ...
	return true;
}

bool CSSimpleSolver::L1Solver(vector<vector<double>> &A, vector<double> &b)
{
	return false;
}

bool CSSimpleSolver::OMPSolver(vector<vector<double>> &A, vector<double> &b)
{
	return false;
}

/*
This solver handles canonical form of LP optimization problem, as follows:
	min CT*x
	s.t. Ax >= b
		  x >= 0;
*/
bool CSSimpleSolver::LPSolver(vector<vector<double>> &A, vector<double> &b, vector<double> &CT, vector<double> &x, double &optimal)
{
	/*
		Input: Constraint matrix A, vector b and objective function CT.
		output: optimal solution x.
		if return = TRUE, x is optimal; otherwise, x may not have meaningful value.
	*/
	if(!A.empty()) {
		// Init
		const int m = A.size(), n = A[0].size();
		lprec *lp;
		int *colno = NULL, i, j, ret = 0;
		REAL *row = NULL;
		lp = make_lp(0, n);		// We will add the rows one by one

		if(lp != NULL) {
			// Move on to load the canonical model
			if(ret == 0) {
				/* create space large enough for one row */
				colno = (int *) malloc(n * sizeof(*colno));
				row = (REAL *) malloc(n * sizeof(*row));
				if((colno == NULL) || (row == NULL))
					ret = 1;
			}
			// Set constraints
			if(ret == 0) {
				set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */
				// Construct the rows of A
				for(i=0; i<m; i++) {
					for(j=0; j<n; j++) {
						colno[j] = j+1;		// May not be efficient in sparse matrix
						row[j] = A[i][j];
					}
					/* add the row to lpsolve */
					if(!add_constraintex(lp, n, row, colno, GE, b[i])) {
						ret = 3;
						printf("Error: Cannot add constrait in %d row of A!", j);
						break;
					}
				}
			}
			// Set objective function
			if(ret == 0) {
				// Construct the objective function 
				set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */
				for(j=0; j<n; j++) {
					colno[j] = j+1;		// May not be efficient in sparse matrix
					row[j] = CT[j];
				}
				 /* set the objective in lpsolve */
				if(!set_obj_fnex(lp, n, row, colno))
				  ret = 4;
			}
			// Call the solver
			if(ret == 0) {
				/* set the object direction to maximize */
				set_minim(lp);

				/* just out of curioucity, now show the model in lp format on screen */
				/* this only works if this is a console application. If not, use write_lp and a filename */
				// write_LP(lp, stdout);		Problematic

				/* I only want to see important messages on screen while solving */
				set_verbose(lp, IMPORTANT);

				/* Now let lpsolve calculate a solution */
				ret = solve(lp);
				if(ret == OPTIMAL)
				  ret = 0;
				else
				  ret = 1;
			}
			// Output optimal solution
			if(ret == 0) {
				/* a solution is calculated, now lets get some results */

				/* objective value */
				//printf("Objective value: %f\n", get_objective(lp));
				optimal = get_objective(lp);

				/* variable values */
				get_variables(lp, row);
				for(j = 0; j < n; j++) x.push_back(row[j]);
				  //printf("%s: %f\n", get_col_name(lp, j + 1), row[j]);

				// Save the model for reuse
				write_lp(lp, "model.lp"); 
				/* we are done now */
			}
			// Clear local memory
			if(row != NULL)
				free(row);
			if(colno != NULL)
				free(colno);
			if(lp != NULL) {
				/* clean up such that all used memory by lpsolve is freed */
				delete_lp(lp);
			}
			return true;
		} else {
			ret = 1;
			printf("Error: Failed to create lp!");
			return false;
		}
	} else {
		std::cout << "Error: Constraint matrix A is empty when solving LP!" << std::endl;
		return false;
	}
}

/*This function solves for Compressed Sensing Model 2, which is Ax + e = b.
	A is m by n matrix, where m>n;
	e is the error vector to be minimized, which is sparse;
	b is the measurement vector.
*/
bool CSSimpleSolver::Solve(vector<vector<double>> &A, vector<double> &b, vector<double> &x, vector<double> &e) 
{
	//// First we need to eliminate A
	//vector<vector<double>> F = leftNullSpace(A);
	//if(F.size() != 0) {
	//	// Then we do the l1-minimization for e, s.t. Fe = Fb.
	//	vector<vector<double>> Y;
	//	vector<vector<double>> B;
	//	B.push_back(b);
	//	multiply_matrix(F, B, Y);

	//	// Then we solve for min|e|, s.t. Fe = Y.
	//	  
	//	

	//} else {
	//	printf("The left nullspace of A is empty!!!\n");
	//	return false;
	//} 

	/* Alternative
		To solve P: min|e|, s.t. Fe = Y;
		We can solve the equivalent problem P': min|e - Ah|;
		Which can be re-expressed as an LP: min 1't, -t <= b-Ag <= t.
		Which can be further denoted as: min [0 1'][g t]'   s.t.   [A I][g t]' >= y && [-A I][g t]' >= -y;
		Let S = [A I; -A I], v = [g t], Y = [y -y]. CT = [0 1'];
	*/

	if(!A.empty()) {
		int m = A.size();
		int n = A[0].size();
		// Construct S
		vector<vector<double>> S;
		for(int i=0; i<m; i++) {
			S.push_back(vector<double>());
			for(int j=0; j<n; j++) {
				S[i].push_back(A[i][j]);
			}
			for(int k=0; k<m; k++) {
				S[i].push_back(i==k);
			}
		}
		for(int i=m; i<2*m; i++) {
			S.push_back(vector<double>());
			for(int j=0; j<n; j++) {
				S[i].push_back(-1*A[i-m][j]);
			}
			for(int k=0; k<m; k++) {
				S[i].push_back((i-m)==k);
			}
		}

		// variable vector v
		vector<double> v(n+m, 0.0);

		// constant vector Y
		vector<double> Y(2*m, 0.0);
		for(int i=0; i<m; i++) {
			Y[i] = b[i];
		}
		for(int i=m; i<2*m; i++) {
			Y[i] = -1*b[i-m];
		}
		
		// Construct objective function
		vector<double> CT(m+n, 0.0);
		for(int i=n; i<m+n; i++) {
			CT[i] = 1.0;
		}

		/*** Solve: min CT*v, s.t. S*v >= Y ***/
		vector<double> solution(m+n, 0.0);
		double optimal;
		bool success = LPSolver(S, Y, CT, solution, optimal);
		if(success) {
			for(int i=0; i<n; i++) {
				x[i] = solution[i];
			}
			for(int i=n; i<m+n; i++) {
				e[i-n] = solution[i];
			}
			return true;
		} else {
			std::cout << "LP solver failed!!" << std::endl;
			return false;
		}
	} else {
		std::cout << "ERROR: A is empty when solving LP!!" << std::endl;
		return false;
	}
}

// Compute Ax = 0, return basis vectors of solution space 
// A is assumed to be sparse matrix
vector<vector<double>> CSSimpleSolver::nullSpace(const vector<vector<double>> &A)
{
	vector<vector<double>> tem;
	if(A.size() != 0) {
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
				std::cout << "Error: Solving nullspace failed!" << std::endl;
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
			//std::cout << "Matrix format error!" << std::endl;
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

		/*print_matrix("U: ", m, n, u, m);
		print_matrix("Singular Values: ", 1, n, s, 1);
		print_matrix("VT: ", n, n, vt, n);*/

		// Define a very small value eps, and set the lower bound to eps*min(m, n), which is similar to Matlab implementation
		//double eps = std::numeric_limits<float>::denorm_min()*max(l_m, l_n);
		double eps = 0.00000000001;
		//printf("eps is %1.11f\n", eps);
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
		//printf("startIndex %d\n", startIndex);
		if(startIndex == n) {
			std::cout << "The nullspace of matrix is empty!" << std::endl;
			return tem;
		}

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
		std::cout << "Null pointer of A in computing nullspace!" << std::endl;
	}
	//std::cout << "tem size: " << tem.size() << std::endl;
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

	if(x.size() != 0) {
		int xm = x.size();
		int xn = x[0].size();
		for(int i=0; i<xn; i++) {
			F.push_back(vector<double>());
			for(int j=0; j<xm; j++) {
				F[i].push_back(x[j][i]);
			}
		}
	} else {
		std::cout << "Error: Compute left null space failed!" << std::endl;
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

bool CSSimpleSolver::multiply_matrix(vector<vector<double>> &A, vector<vector<double>> &B, vector<vector<double>> &C)
{
	if(A.size() == 0 || B.size() == 0) return false;
	int k = A[0].size();
	if(k != B.size()) {
		printf("Cannot multiply by two dismatching matrices,");
		return false;
	} else {
		int m = A.size();
		int n = B[0].size();
		for(int mm = 0; mm<m; mm++) {
			C.push_back(vector<double>());
			for(int nn=0; nn<n; nn++) {
				double sum = 0.0;
				for(int kk=0; kk<k; kk++) {
					sum += A[mm][kk]*B[kk][nn];
				}
				C[mm].push_back(sum); sum = 0.0;
			}
		}
		return true;
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
	/*vector<vector<double>> A;
	double a0[] = {-1, 2};
	double a1[] = {1, 0};
	double a2[] = {2, 1};
	double a3[] = {4, -7};
	vector<double> A0(a0, a0 + sizeof(a0)/sizeof(double));
	vector<double> A1(a1, a1 + sizeof(a1)/sizeof(double));
	vector<double> A2(a2, a2 + sizeof(a2)/sizeof(double));
	vector<double> A3(a3, a3 + sizeof(a3)/sizeof(double));
	A.push_back(A0); A.push_back(A1); A.push_back(A2); A.push_back(A3);*/


	CSSimpleSolver css;
	/*vector<vector<double>> x = css.leftNullSpace(A);
	for(int i=0; i<x.size(); i++) {
		for(int j=0; j<x[0].size()-1; j++) {
			std::cout << x[i][j] << "\t" ;
		}
		std::cout << x[i][x[0].size()-1] << std::endl;
	}*/

	// Test matrix multiplication
	/*vector<vector<double>> B;
	double b0[] = {1};
	double b1[] = {2};
	vector<double> B0(b0, b0 + sizeof(b0)/sizeof(double));
	vector<double> B1(b1, b1 + sizeof(b1)/sizeof(double));
	B.push_back(B0); B.push_back(B1);
	vector<vector<double>> C;
	css.multiply_matrix(A, B, C);
	for(int i=0; i<C.size(); i++) {
		for(int j=0; j<C[0].size(); j++) {
			std::cout << C[i][j] << "\t";
		}
		std::cout << std::endl;
	}*/

	// Test LPSolver()
	double a0[] = {-120, -210};
	double a1[] = {-110, -30};
	double a2[] = {-1, -1};
	vector<vector<double>> A;
	vector<double> A0(a0, a0 + sizeof(a0)/sizeof(double));
	vector<double> A1(a1, a1 + sizeof(a1)/sizeof(double));
	vector<double> A2(a2, a2 + sizeof(a2)/sizeof(double));
	A.push_back(A0); A.push_back(A1); A.push_back(A2);

	double b[] = {-15000, -4000, -75};
	vector<double> B(b, b + sizeof(b)/sizeof(double));

	double ct[] = {-143, -60};
	vector<double> CT(ct, ct + sizeof(ct)/sizeof(double));

	vector<double> x;
	double optimal = 0.0;
	css.LPSolver(A, B, CT, x, optimal);
	for(int i=0; i<x.size(); i++) printf("%f\t", x[i]);
	printf("\nOptimal Solution: %f\n", optimal);

	return 0;
}
	