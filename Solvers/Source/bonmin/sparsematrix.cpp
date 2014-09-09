// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#include "sparsematrix.hpp"
#include "matlabexception.hpp"
#include "iterate.hpp"

// Function definitions for class SparseMatrix.
// ---------------------------------------------------------------
SparseMatrix::SparseMatrix (const mxArray* ptr) 
  : jc(0), ir(0), x(0), linA(0), vec_c(0) {

    // Get the height, width and number of non-zeros.
    h   = (int) mxGetM(ptr);
    w   = (int) mxGetN(ptr);
    nnz = getSizeOfSparseMatrix(ptr);  
    
    // Deep Copy the row and column indices, and the values of the nonzero entries.
    jc = new mwIndex[w+1];
    ir = new mwIndex[nnz];
    x  = new double[nnz];
    copymemory(mxGetJc(ptr),jc,w+1); 
    copymemory(mxGetIr(ptr),ir,nnz); 
    copymemory(mxGetPr(ptr),x,nnz);  
    // Copy the pointer to the MATLAB matrix (for Matrix*Vector later (if required))
    linA = (mxArray*)ptr;
    //Set Matlab x memory to null
    prhs[1] = NULL;
}
  
 //Copy constructor
 SparseMatrix::SparseMatrix (const SparseMatrix *obj)
  : jc(0), ir(0), x(0), linA(0), vec_c(0) {      
    // Copy height, width and number of non-zeros.
    h   = obj->h;
    w   = obj->w;
    nnz = obj->nnz;

    // Deep Copy the row and column indices, and the values of the nonzero entries.
    jc = new mwIndex[w+1];
    ir = new mwIndex[nnz];
    x  = new double[nnz];
    copymemory(obj->jc,jc,w+1); 
    copymemory(obj->ir,ir,nnz); 
    copymemory(obj->x,x,nnz); 
    //Set Matlab x memory to null
    prhs[1] = NULL;
}
 
//Preallocation Constructor
SparseMatrix::SparseMatrix (int h_, int w_, int nnz_) 
 : jc(0), ir(0), x(0), linA(0), vec_c(0) {     
    // Copy height, width and number of non-zeros.
    this->h   = h_;
    this->w   = w_;
    this->nnz = nnz_;

    //Create Index and Element Memory
    jc = new mwIndex[w+1];
    ir = new mwIndex[nnz];
    x  = new double[nnz]; 
    //Set Matlab x memory to null
    prhs[1] = NULL;
} 

SparseMatrix::~SparseMatrix() {
    if (jc)     {delete[] jc;   jc=NULL;}
    if (ir)     {delete[] ir;   ir=NULL;}
    if (x)      {delete[] x;    x=NULL;}
    if(vec_c)   {mxDestroyArray(vec_c);   vec_c=NULL;}
    if(prhs[1]) {mxDestroyArray(prhs[1]); prhs[1]=NULL;}
}

size_t SparseMatrix::numelems (int c) const {
  return jc[c+1] - jc[c];
}

void SparseMatrix::getColsAndRows (int* cols, int* rows) const {
       
  // Repeat for each column in the matrix, then repeat for each
  // non-zero entry in the current column.
  for (int c = 0, i = 0; c < w; c++)
    for ( ; i < (int) jc[c+1]; i++) {
      cols[i] = (int) c;
      rows[i] = (int) ir[i];
    }
}

bool SparseMatrix::copyto (SparseMatrix& dest) const {
  bool match, matchrow, matchcol;  // Loop variables.

  // Initialize the destination values to zero, because we might not
  // have corresponding nonzero entries from the source for some of
  // the destination nonzero entries.
  for (int i = 0; i < (int)dest.nnz; i++)
    dest.x[i] = 0;

  // In order to properly copy the elements from one sparse matrix to
  // another, the non-zero elements of the destination matrix must be
  // a superset of the collection of non-zero elements in the source
  // matrix. If not, the copy operation is invalid, and the function
  // returns false. Repeat for each column.
  int i = 0;  // Index of element in source.
  int j = 0;  // Index of element in destination.
  for (int c = 0; c < dest.w; c++) {

    // Repeat for each non-zero element in the destination column.
    for ( ; j < (int) dest.jc[c+1]; j++) {

      // Let's check to see if the source column matches the
      // destination column, and the source row matches the
      // destination row. The first line checks to see if the row
      // indices match. The second line checks to see if the column
      // indices match. (If the source matrix is valid, there is no
      // need to check whether the column index of the source entry is
      // LESS THAN the destination column index since we assume that
      // there are less entries in the source, and under that
      // assumption we move faster through the source matrix than we
      // do through the destination matrix.)
      matchrow = (ir[i] == dest.ir[j]);
      matchcol = (i >= (int) jc[c]) && (i < (int) jc[c+1]);
      match    = matchrow && matchcol;

      // If the row & column indices match, then we copy the source
      // entry value and move on to the next non-zero entry in the
      // source.
      dest.x[j] = match * x[i];
      i        += match;
    }      
  }

  // If we've reached the end of the loop and we haven't visited all
  // the non-zero entries in the source matrix, it means that the
  // source matrix is invalid, and we throw an error.
  if (i < (int)nnz)
    return false;
  else
    return true;
}

void SparseMatrix::copyto (double* dest) const {
  copymemory(x,dest,nnz);
}

//Sparse Matrix*Vector using MATLAB (assumes constant structure as per linear A)
void SparseMatrix::SpMatrixVec(const Iterate& xin, double *c) {
    //Sizes
    int M = this->h;
    int N = this->w;
    //Check dims    
    if(w != numvars(xin))
        throw MatlabException("To multiply a sparse matrix by a vector the number of columns in the matrix must equal the number of rows in the vector");
    
    if(linA == NULL)
        throw MatlabException("Error with sparse matrix linear A memory");
    
    //Check if existing vec_c memory exists, if not, create it and assign pointers at the same time
    if(vec_c == NULL) {    
        prhs[0] = this->linA;
        prhs[1] = mxCreateDoubleMatrix(N,1,mxREAL);
        vec_c = mxCreateDoubleMatrix(M,1,mxREAL);
    }
    //Copy in current x iterate to MATLAB memory
    xin.copyto(mxGetPr(prhs[1]));
    //Call Matlab to evaluate sparse Matrix * Vector
    try {
        mexCallMATLAB(1,&vec_c,2,prhs,"mtimes");
    }
    catch (std::exception ME) {
        const char* what = ME.what();
        if (what) {
            mexPrintf("Matlab exception:\n%s", what);      
            throw MatlabException("There was an error when executing the linear constraints");
        }
    }
    catch (...) {
        throw MatlabException("There was an error when executing the linear constraints");
    }
    //Copy results to output c pointer
    memcpy(c,mxGetPr(vec_c),M*sizeof(double));
}

//Vertical Concatentation of two sparse matrices into a new matrix
void SparseMatrix::VertConcatenate(const SparseMatrix *obj, SparseMatrix *SpCat) {
    size_t ind = 0, i = 0, j = 0, k = 0;
    //Check dims
    if(this->w != obj->w)
        throw MatlabException("To vertically concatenate sparse matrices both must have the same number of columns");
    
    //Save Sizes
    int N = this->w;
    int M1 = this->h;
    int M2 = obj->h;
    
    //Check sizes
    if(height(*SpCat) != M1+M2)
        throw MatlabException("Wrong number of rows in resulting matrix after concatenation");
    if(width(*SpCat) != N)
        throw MatlabException("Wrong number of columns in resulting matrix after concatenation");
    if(SpCat->numelems() != this->numelems()+obj->numelems())
        throw MatlabException("Wrong number of nnzs in resulting matrix after concatenation");
    
    //Concatenate into new matrix
    SpCat->jc[0] = 0;
    for(i = 1; i <= (size_t)N; i++) {
        SpCat->jc[i] = this->jc[i] + obj->jc[i];
        for(j = 0; j < this->jc[i]-this->jc[i-1]; j++) {
            ind = this->jc[i-1]+j;
            SpCat->ir[k] = this->ir[ind];
            SpCat->x[k++] = this->x[ind];
        }
        for(j = 0; j < obj->jc[i]-obj->jc[i-1]; j++) {
            ind = obj->jc[i-1]+j;
            SpCat->ir[k] = obj->ir[ind]+M1;
            SpCat->x[k++] = obj->x[ind];
        }
    }
}

// Function definitions for static members of class SparseMatrix.
// -----------------------------------------------------------------
size_t SparseMatrix::getSizeOfSparseMatrix (const mxArray* ptr) {
  
  // Get the width (the number of columns) of the matrix.
  size_t w = (size_t) mxGetN(ptr);
  
  // The precise number of non-zero elements is contained in the
  // last entry of the jc array. (There is one jc entry for each
  // column in the matrix, plus an extra one.)
  mwIndex* jc = mxGetJc(ptr);
  return jc[w];    
}  

bool SparseMatrix::isLowerTri (const mxArray* ptr) {

  // Get the height and width of the matrix.
  int h = (int) mxGetM(ptr);
  int w = (int) mxGetN(ptr);
  
  // Check whether the sparse matrix is symmetric.
  if (h != w)
    return false;

  // A sparse lower triangular matrix has the property that the
  // column indices are always less than or equal to the row
  // indices.
  bool     b  = true;  // The return value.
  mwIndex* jc = mxGetJc(ptr);
  mwIndex* ir = mxGetIr(ptr);
  for (int c = 0, i = 0; c < w; c++)
    for ( ; i < (int) jc[c+1]; i++)
      b = b && (c <= (int) ir[i]);
  
  return b;
}

bool SparseMatrix::inIncOrder (const mxArray* ptr) {
  bool     b  = true;
  int      w  = (int) mxGetN(ptr);
  mwIndex* jc = mxGetJc(ptr);
  mwIndex* ir = mxGetIr(ptr);

  for (int c = 0, i = 0; c < w; c++)
    if (jc[c+1] > jc[c]) {
      i++;
      for ( ; i < (int) jc[c+1]; i++)
	b = b && (ir[i] > ir[i-1]);
    }

  return b;
}
