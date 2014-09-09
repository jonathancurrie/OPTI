#include "lpsolvecaller.h"

#if defined STANDALONE
#pragma comment(lib, "libmx.lib")
#pragma comment(lib, "libmat.lib")
#pragma comment(lib, "libmex.lib")
#pragma comment(lib, "libmx.lib")
#pragma comment(lib, "libeng.lib")
/* #pragma comment(lib, "libmatlb.lib") */
#endif

//JONNY EDIT - doesn't work for 64bit!
// #ifndef SINCE73 /* NEW MATLAB 7.3 SIZE POINTERS */
// # define mwSize int
// # define mwIndex int
// #endif

/* Function to get a real scalar with error checking */

Double GetRealScalar(structlpsolvecaller *lpsolvecaller, int element)
{
        rMatrix pm = GetpMatrix(lpsolvecaller, element);

	if ((GetM(lpsolvecaller, pm) == 1) && (GetN(lpsolvecaller, pm) == 1)
	  && (IsNumeric(pm)) && (!IsComplex(pm)) ) {
		return(mxGetScalar(pm));
	} else {
		ErrMsgTxt(lpsolvecaller, "Expecting a scalar argument.");
	}
        return(0.0);
}

#define GetVector(lpsolvecaller, element, vec, cast, start, len, exactcount, ret) \
{ \
	int	j, k, k1, count = 0; \
        mwSize  m, n; \
        mwIndex	*ir, *jc; \
	double	*pr, *pr0; \
        cast    *vec0; \
        rMatrix pm = GetpMatrix(lpsolvecaller, element); \
 \
	m = GetM(lpsolvecaller, pm); \
	n = GetN(lpsolvecaller, pm); \
 \
	if ( !((m == 1) || (n == 1)) || \
             ((m == 1) && (((exactcount) && (len != n)) || ((!exactcount) && (n > len)))) || \
             ((n == 1) && (((exactcount) && (len != m)) || ((!exactcount) && (m > len)))) || \
	     !IsNumeric(pm) || IsComplex(pm)) { \
		ErrMsgTxt(lpsolvecaller, "invalid vector."); \
	} \
 \
	pr = GetPr(pm); \
 \
	if (!IsSparse(pm)) { \
                if (n == 1) \
                        len = m; \
                else \
                        len = n; \
                vec += start; \
		for (k = 0; k < len; k++, pr++, vec++) { \
			*vec = (cast) *pr; \
		} \
                count = len; \
	} else if (IsSparse(pm)) { \
		jc = mxGetJc(pm); \
		ir = mxGetIr(pm); \
                pr0 = pr; \
                vec0 = vec; \
                for (j = 0; j < n; j++) { \
                        k = jc[j]; \
                        k1 = jc[j + 1]; \
                        pr = pr0 + k; \
                        vec = vec0; \
                        vec += start + j * m; \
                        for (; k < k1; k++, pr++) { \
				vec[ir[k]] = (cast) *pr; \
                                count++; \
			} \
		} \
	} else { \
		ErrMsgTxt(lpsolvecaller, "Can't figure out this matrix."); \
	} \
 \
        ret = count; \
}

/* Functions to get len elements from a MATLAB vector. Matrix
   can be either full or sparse. Elements are stored in indices
   start..start+n-1  Errors out if the MATLAB vector is not length len */

int GetIntVector(structlpsolvecaller *lpsolvecaller, int element, int *vec, int start, int len, int exactcount)
{
        int ret;

        GetVector(lpsolvecaller, element, vec, int, start, len, exactcount, ret);

        return(ret);
}

int GetRealVector(structlpsolvecaller *lpsolvecaller, int element, Double *vec, int start, int len, int exactcount)
{
        int ret;

        GetVector(lpsolvecaller, element, vec, Double, start, len, exactcount, ret);

        return(ret);
}


/* Function to get max len elements from a MATLAB sparse vector. Matrix
   can be either full or sparse. Elements are stored in indices
   start..start+n-1  Errors out if the MATLAB vector is longer than length len */

int GetRealSparseVector(structlpsolvecaller *lpsolvecaller, int element, Double *vec, int *index, int start, int len, int col)
{
	int	j, k, k1, start1, count = 0;
        mwSize  m, n;
	mwIndex	*ir, *jc;
	double	*pr, *pr0;
        rMatrix pm = GetpMatrix(lpsolvecaller, element);

	m = GetM(lpsolvecaller, pm);
	n = GetN(lpsolvecaller, pm);

	if (  ((col == 0) && (((m != 1) && (n != 1)) || ((m == 1) && (n > len)) || ((n == 1) && (m > len)))) ||
              ((col != 0) && ((m > len) || (col > n))) ||
	      !IsNumeric(pm) ||
              IsComplex(pm)  ) {
		ErrMsgTxt(lpsolvecaller, "invalid vector.");
	}

	pr = GetPr(pm);

	if (!IsSparse(pm)) {
                if ((((n == 1) || (col != 0)) && (m != len)) || ((col == 0) && (m == 1) && (n != len)))
                        ErrMsgTxt(lpsolvecaller, "invalid vector.");

                if (col)
                	pr += (col - 1) * m;
                for (k = 0; k < len; k++, pr++) {
                        if (*pr) {
				*(vec++) = *pr;
                        	*(index++) = start + k;
                        	count++;
                        }
		}
	} else if (IsSparse(pm)) {
                int j1, j2;

		jc = mxGetJc(pm);
		ir = mxGetIr(pm);
                pr0 = pr;
                if (col == 0) {
                        j1 = 0;
                        j2 = n;
                }
                else {
                        j1 = col - 1;
                        j2 = col;
                }
		for (j = j1; j < j2; j++) {
                        k = jc[j];
                        k1 = jc[j + 1];
                        pr = pr0 + k;
                        start1 = start;
                        if (col == 0)
                        	start1 += j * m;
                        for (; k < k1; k++, pr++, vec++, index++) {
                                *vec = *pr;
                                *index = start1 + (int)ir[k]; //JONNY EDIT
                                count++;
			}
		}
	} else {
		ErrMsgTxt(lpsolvecaller, "Can't figure out this matrix.");
	}

        return(count);
}


int GetString(structlpsolvecaller *lpsolvecaller, pMatrix ppm, int element, char *buf, int size, int ShowError)
{
        rMatrix pm;

        if (ppm != NULL)
                pm = ppm[element];
        else
                pm = GetpMatrix(lpsolvecaller, element);

	if (!mxIsChar(pm)) {
                if (ShowError)
                	ErrMsgTxt(lpsolvecaller, "Expecting a character element.");
                return(FALSE);
        }
	mxGetString(pm, buf, size);
        return(TRUE);
}

strArray GetCellCharItems(structlpsolvecaller *lpsolvecaller, int element, int len, int ShowError)
{
        int  i;
        mwSize m, n;
        rMatrix pm = GetpMatrix(lpsolvecaller, element);
        pMatrix pa0 = NULL, **pa;

        if (!mxIsCell(pm)) {
                if (ShowError)
                        ErrMsgTxt(lpsolvecaller, "Expecting a cell argument.");
        }
        else {
                m = GetM(lpsolvecaller, pm);
                n = GetN(lpsolvecaller, pm);
                if (!(((m == 1) && (n == len)) || ((n == 1) && (m == len))))
                        ErrMsgTxt(lpsolvecaller, "invalid vector.");
                pa = pa0 = (pMatrix) matCalloc(len, sizeof(*pa));
                for (i = 0; i < len; i++) {
                	*pa = mxGetCell(pm, i);
                	if (!mxIsChar(*pa))
                                break;
                        pa++;
                }
                if (i < len) {
                        matFree(pa0);
                	ErrMsgTxt(lpsolvecaller, "Expecting a character cell element.");
                }
        }
        return(pa0);
}


void FreeCellCharItems(strArray pa, int len)
{
	matFree(pa);
}

double *CreateDoubleMatrix(structlpsolvecaller *lpsolvecaller, int m, int n, int element)
{
	lpsolvecaller->plhs[element] = mxCreateDoubleMatrix(m, n, mxREAL);
	return(GetPr(lpsolvecaller->plhs[element]));
}

double *CreateDoubleSparseMatrix(structlpsolvecaller *lpsolvecaller, int m, int n, int element)
{
        int nzmax;

        if (m > n)
                nzmax = m;
        else
                nzmax = n;
        if (nzmax < 10)
                nzmax = 10;

        lpsolvecaller->plhs[element] = mxCreateSparse(m, n, nzmax, mxREAL);

        return(NULL);
}

void SetColumnDoubleSparseMatrix(structlpsolvecaller *lpsolvecaller, int element, int m, int n, double *mat, int column, double *arry, int *index, int size, int *nz)
{
        int i, ii, nzmax;
        mwIndex  *jcs, *irs;
        double *sr, a;

        jcs = mxGetJc(lpsolvecaller->plhs[element]);
        irs = mxGetIr(lpsolvecaller->plhs[element]);
        sr  = mxGetPr(lpsolvecaller->plhs[element]);
        if (m > n)
                nzmax = m;
        else
                nzmax = n;
        if (nzmax < 10)
                nzmax = 10;

        jcs[column - 1] = *nz;
        for (i = 0; (i < size); i++) {
                a = arry[i];
                if (a) {
                        if (index == NULL)
                                ii = i;
                        else
                                ii = index[i] - 1;
                        if ((*nz != 0) && (((*nz) % nzmax) == 0)) {
                                mxSetNzmax(lpsolvecaller->plhs[element], *nz + nzmax);
                                mxSetPr(lpsolvecaller->plhs[element], matRealloc(sr, (*nz + nzmax) * sizeof(double)));
                                mxSetIr(lpsolvecaller->plhs[element], matRealloc(irs, (*nz + nzmax) * sizeof(mwIndex)));
                                jcs = mxGetJc(lpsolvecaller->plhs[element]);
                                irs = mxGetIr(lpsolvecaller->plhs[element]);
                                sr  = mxGetPr(lpsolvecaller->plhs[element]);
                        }
                        sr[*nz] = a;
                        irs[*nz] = ii;
                        (*nz)++;
                }
        }
        jcs[column] = *nz;
}

void CreateString(structlpsolvecaller *lpsolvecaller, char **str, int n, int element)
{
        rMatrix pa;

        if (n == 1) {
         	pa = mxCreateString(str[0]);
                if (lpsolvecaller->plhs != NULL)
                	lpsolvecaller->plhs[element] = pa;
        }
        else {
                int i;
        	rMatrix pm;

                pm = lpsolvecaller->plhs[element] = mxCreateCellMatrix(1, n);
                for (i = 0; i < n; i++) {
                        pa = mxCreateString(str[i]);
                	mxSetCell(pm, i, pa);
        	}
        }
}
