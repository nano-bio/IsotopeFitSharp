// SPQR.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
//#include <SuiteSparseQR_C.h>
#include <iostream>
#include <sstream>
#include <SuiteSparseQR.hpp>

#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>

extern "C"
{

	__declspec(dllexport) int Pica(
		void ** handle,
		int ord, double tol,
		int rows, int cols, int nzCount, double vals[], __int64 rowInd[], __int64 colPtr[],
		int * RRows, int * RCols, int * RNzCount, double ** RVals, __int64 ** RRowInd, __int64 ** RColPtr, __int64 ** ROrdering)
	{
		return 0;
	}

#ifdef _WIN64
	__declspec(dllexport) int SparseQR(
		void ** handle,
		int ord, double tol,
		int rows, int cols, int nzCount, double vals[], __int64 rowInd[], __int64 colPtr[],
		int * RRows, int * RCols, int * RNzCount, double ** RVals, __int64 ** RRowInd, __int64 ** RColPtr, __int64 ** ROrdering)
#else
	__declspec(dllexport) int SparseQR(
		void ** handle,
		int ord, double tol,
		int rows, int cols, int nzCount, double vals[], __int64 rowInd[], __int64 colPtr[],
		int * RRows, int * RCols, int * RNzCount, double ** RVals, __int64 ** RRowInd, __int64 ** RColPtr, long ** ROrdering)
#endif
	{
				// create workspace objects and add their addresses to the handle array
				cholmod_common *cc = new cholmod_common();
				cholmod_sparse *A = new cholmod_sparse();
				cholmod_sparse *R = new cholmod_sparse();
		
				// this is what we want the handle to look like
				handle[0] = cc;
				handle[1] = A;
				handle[2] = R;
				
				// internal hard-coded parameters
		#ifdef _WIN64
				__int64 econ = cols;			// how many rows of R will be returned - economic mode
		#else
				long econ = cols;
		#endif
				bool sorted = true;				// values in rowInd array are sorted within the columns
				bool packed = true;				// columns are stored continuously, without zeroes
				int matrixType = 0;				// matrix is unsymmetric
		
				// initialize the workspace
				cholmod_l_start(cc);
		
				// create the input matrix object and fill it with supplied values
				A = cholmod_l_allocate_sparse(rows, cols, nzCount, sorted, packed, matrixType, CHOLMOD_REAL, cc);
				A->x = vals;
				A->i = rowInd;
				A->p = colPtr;
				A->nzmax = nzCount;
		
				// check for input validity
				int check = cholmod_l_check_sparse(A, cc);	
				if (check != 1) return -1;
		
				// QR factorization
				__int64 rank = SuiteSparseQR<double>(ord, tol, econ, A, &R, ROrdering, cc);	//TODO: ordering?
		
				// set outputs
				*RVals = (double*)R->x;
				*RRowInd = (__int64*)R->i;
				*RColPtr = (__int64*)R->p;
				*RRows = R->ncol;
				*RCols = R->nrow;
				*RNzCount = R->nzmax;

		return 0;
	}

	__declspec(dllexport) void SparseQRDispose(void ** handle)
	{
		cholmod_common *cc = (cholmod_common*)handle[0];
		cholmod_sparse *A = (cholmod_sparse*)handle[1];
		cholmod_sparse *R = (cholmod_sparse*)handle[2];

		cholmod_l_free_sparse(&A, cc);
		cholmod_l_free_sparse(&R, cc);
		cholmod_l_finish(cc);
		cholmod_l_free_work(cc);

		delete cc;
		delete A;
		delete R;
	}

	//// X = A\dense(B) using default ordering and tolerance
	//template <typename Entry> cholmod_dense *SuiteSparseQR
	//(
	//	cholmod_sparse *A,      // m-by-n sparse matrix
	//	cholmod_dense  *B,      // m-by-nrhs
	//	cholmod_common *cc      // workspace and parameters
	//);

	__declspec(dllexport) int SparseSolve(
		void ** handle,
		int rows, int cols, int nzCount, double * vals, __int64 * rowInd, __int64 * colPtr,
		double * bVals,
		double ** xVals)
	{
		//_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
		//_CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_DEBUG);

		std::stringstream oss;
		std::streambuf* old = std::cerr.rdbuf(oss.rdbuf());

		// create workspace objects and add their addresses to the handle array
		cholmod_common *cc = new cholmod_common();
		cholmod_sparse *A = new cholmod_sparse();
		cholmod_dense *b = new cholmod_dense();
		cholmod_dense *x = new cholmod_dense();

		// this is what we want the handle to look like
		handle[0] = cc;
		handle[1] = A;
		handle[2] = b;
		handle[3] = x;

		// internal hard-coded parameters
		bool sorted = true;				// values in rowInd array are sorted within the columns
		bool packed = true;				// columns are stored continuously, without zeroes
		int matrixType = 0;				// matrix is unsymmetric

		// initialize the workspace
		cholmod_l_start(cc);

		// create the input matrix object and fill it with supplied values
		A = cholmod_l_allocate_sparse(rows, cols, nzCount, sorted, packed, matrixType, CHOLMOD_REAL, cc);
		A->x = vals;
		A->i = rowInd;
		A->p = colPtr;
		//A->nzmax = nzCount;

		// check for input validity
		int check = cholmod_l_check_sparse(A, cc);
		if (check != 1) return -1;

		b = cholmod_l_allocate_dense(rows, 1, rows, CHOLMOD_REAL, cc);
		b->x = bVals;
		
		check = cholmod_l_check_dense(b, cc);
		if (check != 1) return -1;

		// system solution
		x = SuiteSparseQR<double>(A, b, cc);

		// set outputs
		*xVals = (double*)x->x;

		return x->nzmax;
	}

	__declspec(dllexport) void SparseSolveDispose(void ** handle)
	{
		//_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
		//_CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_DEBUG);

		cholmod_common *cc = (cholmod_common*)handle[0];
		cholmod_sparse *A = (cholmod_sparse*)handle[1];
		cholmod_dense *b = (cholmod_dense*)handle[2];
		cholmod_dense *x = (cholmod_dense*)handle[3];

		cholmod_l_free_sparse(&A, cc);
		cholmod_l_free_dense(&b, cc);
		cholmod_l_free_dense(&x, cc);
		cholmod_l_free_work(cc);
		cholmod_l_finish(cc);

		delete cc;
		delete A;
		delete b;
		delete x;
	}
}
