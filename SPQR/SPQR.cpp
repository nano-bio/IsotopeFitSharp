// SPQR.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
//#include <SuiteSparseQR_C.h>
//#include <iostream>
#include <SuiteSparseQR.hpp>

extern "C"
{
	__declspec(dllexport) int SparseQR(
		void ** handle,
		int ord, double tol,
		int rows, int cols, int nzCount, double vals[], __int64 rowInd[], __int64 colPtr[],
		int * RRows, int * RCols, int * RNzCount, double ** RVals, __int64 ** RRowInd, __int64 ** RColPtr, __int64 ** ROrdering)
	{
		// create workspace objects and add their addresses to the handle array
		cholmod_common *cc = new cholmod_common();
		cholmod_sparse *A = new cholmod_sparse();
		cholmod_sparse *R = new cholmod_sparse();

		handle[0] = cc;
		handle[1] = A;
		handle[2] = R;
		
		// internal hard-coded parameters
		__int64 econ = cols;			// how many rows of R will be returned - economic mode
		bool sorted = true;				// values in rowInd array are sorted per column
		bool packed = true;				// columns are stored continuously
		int matrixType = 0;				// matrix is unsymmetric

		// initialize the workspace
		cholmod_l_start(cc);

		// create the input matrix object and fill it with supplied values
		A = cholmod_l_allocate_sparse(rows, cols, nzCount, sorted, packed, matrixType, CHOLMOD_REAL, cc);

		A->x = vals;
		A->i = rowInd;
		A->p = colPtr;
		A->nzmax = nzCount;

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
}
