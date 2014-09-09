// Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         September 25, 2008

#include "options.hpp"
#include "matlabexception.hpp"

// Function definitions for class Options.
// -----------------------------------------------------------------
Options::Options (const Iterate& x, Ipopt::IpoptApplication& app, 
		  const mxArray* ptr) 
  : n(numvars(x)), m(0), nlin(0), nnlin(0), lb(0), ub(0), cl(0), cu(0), zl(0), zu(0),
    lambda(0), auxdata(0),

    // Process the IPOPT options.
    ipopt(app,mxGetField(ptr,0,"ipopt")) { 

  double neginfty = ipopt.getNegInfty();  // Negative infinity.
  double posinfty = ipopt.getPosInfty();  // Positive infinity.

  // Load the bounds on the variables.
  lb = loadLowerBounds(n,ptr,neginfty);
  ub = loadUpperBounds(n,ptr,posinfty);

  // Load the bounds on the constraints.
  m = loadConstraintBounds(ptr,cl,cu,neginfty,posinfty,nlin,nnlin);

  // Load the Lagrange multipliers.
  loadMultipliers(n,m,ptr,zl,zu,lambda);

  // Load the auxiliary data.
  auxdata = mxGetField(ptr,0,"auxdata");
}

Options::~Options() {
  if (lb) delete[] lb;
  if (ub) delete[] ub;
  if (cl) delete[] cl;
  if (cu) delete[] cu;
}

// Function definitions for static members of class Options.
// -----------------------------------------------------------------
double* Options::loadLowerBounds(int n, const mxArray* ptr, double neginfty) {
  double*        lb;  // The return value.
  const mxArray* p  = mxGetField(ptr,0,"lb");

  if (p) {

    // Load the upper bounds and check to make sure they are valid.
    int N = Iterate::getMatlabData(p,lb);
    if (N != n)
      throw MatlabException("Lower bounds array must have one element for \
each optimization variable");

    // Convert MATLAB's convention of infinity to IPOPT's convention
    // of infinity.
    for (int i = 0; i < n; i++)
      if (mxIsInf(lb[i]))
	lb[i] = neginfty;
  } else {
    
    // If the lower bounds have not been specified, set them to
    // negative infinity.
    lb = new double[n];
    for (int i = 0; i < n; i++)
      lb[i] = neginfty;
  }

  return lb;
}

double* Options::loadUpperBounds(int n, const mxArray* ptr, double posinfty) {
  double* ub;  // The return value.

  // Load the upper bounds on the variables.
  const mxArray* p = mxGetField(ptr,0,"ub");
  if (p) {

    // Load the upper bounds and check to make sure they are valid.
    int N = Iterate::getMatlabData(p,ub);
    if (N != n)
      throw MatlabException("Upper bounds array must have one element for \
each optimization variable");

    // Convert MATLAB's convention of infinity to IPOPT's convention
    // of infinity.
    for (int i = 0; i < n; i++)
      if (mxIsInf(ub[i]))
	ub[i] = posinfty;
  } else {

    // If the upper bounds have not been specified, set them to
    // positive infinity.
    ub = new double[n];
    for (int i = 0; i < n; i++)
      ub[i] = posinfty;
  }

  return ub;
}

int Options::loadConstraintBounds (const mxArray* ptr, double*& cl, 
				   double*& cu, double neginfty,
				   double posinfty, int &lin, int &nlin) {
  int m = 0;  // The return value is the number of constraints.
  int tm;

  // LOAD CONSTRAINT BOUNDS.
  // If the user has specified constraints bounds, then she must
  // specify *both* the lower and upper bounds.
  const mxArray* pl = mxGetField(ptr,0,"cl");
  const mxArray* pu = mxGetField(ptr,0,"cu");
  const mxArray* prl = mxGetField(ptr,0,"rl"); //linear constraint bounds
  const mxArray* pru = mxGetField(ptr,0,"ru");
  if (pl || pu || prl || pru) {

    // Check to make sure the constraint bounds are valid.
    if ((!pl ^ !pu) || (!prl ^ !pru))
      throw MatlabException("You must specify both lower and upper bounds on the constraints");
    if (pl && (!mxIsDouble(pl) || !mxIsDouble(pu) || (mxGetNumberOfElements(pl) != mxGetNumberOfElements(pu))))
      throw MatlabException("The nonlinear constraints lower and upper bounds must both be double-precision arrays with the same number of elements");
    if (prl && (!mxIsDouble(prl) || !mxIsDouble(pru) || (mxGetNumberOfElements(prl) != mxGetNumberOfElements(pru))))
      throw MatlabException("The linear constraints lower and upper bounds must both be double-precision arrays with the same number of elements");
    // Get the number of constraints.
    if(pl && prl) {
        lin = (int)mxGetNumberOfElements(prl);
        nlin = (int)mxGetNumberOfElements(pl);
        m = lin+nlin;
    }
    else if(pl) {
        lin = 0;
        nlin = (int)mxGetNumberOfElements(pl);
        m = nlin;        
    }
    else {
        lin = (int)mxGetNumberOfElements(prl);
        nlin = 0;
        m = lin;
    }

    // Load the lower bounds on the constraints and convert MATLAB's
    // convention of infinity to IPOPT's convention of infinity.
    cl = new double[m];
    cu = new double[m];
    if(pl && prl) {
        tm = (int)mxGetNumberOfElements(pl);
        copymemory(mxGetPr(pl),cl,tm);
        copymemory(mxGetPr(pu),cu,tm);
        copymemory(mxGetPr(prl),&cl[tm],(int)mxGetNumberOfElements(prl));
        copymemory(mxGetPr(pru),&cu[tm],(int)mxGetNumberOfElements(pru));
    }
    else if(pl) {
        copymemory(mxGetPr(pl),cl,m);
        copymemory(mxGetPr(pu),cu,m);
    }
    else {
        copymemory(mxGetPr(prl),cl,m);
        copymemory(mxGetPr(pru),cu,m);
    }

    // Convert MATLAB's convention of infinity to IPOPT's convention
    // of infinity.
    for (int i = 0; i < m; i++) {
      if (mxIsInf(cl[i])) cl[i] = neginfty;
      if (mxIsInf(cu[i])) cu[i] = posinfty;
    }
  }

  return m;
}

void Options::loadMultipliers (int n, int m, const mxArray* ptr, 
			       double*& zl, double*& zu, double*& lambda) {
  const mxArray* p;
  
  // Load the Lagrange multipliers associated with the lower bounds.
  p = mxGetField(ptr,0,"zl");
  if (p) {
    if (!mxIsDouble(p) || (int) mxGetNumberOfElements(p) != n)
      throw MatlabException("The initial point for the Lagrange multipliers \
associated with the lower bounds must be a double-precision array with one \
element for each optimization variable");
    zl = mxGetPr(p);
  } else
    zl = 0;

  // Load the Lagrange multipliers associated with the upper bounds.
  p = mxGetField(ptr,0,"zu");
  if (p) {
    if (!mxIsDouble(p) || (int) mxGetNumberOfElements(p) != n)
      throw MatlabException("The initial point for the Lagrange multipliers \
associated with the upper bounds must be a double-precision array with one \
element for each optimization variable");
    zu = mxGetPr(p);
  } else
    zu = 0;

  // Load the Lagrange multipliers associated with the equality and
  // inequality constraints.
  p = mxGetField(ptr,0,"lambda");
  if (p) {
    if (m>0 && (!mxIsDouble(p) || (int) mxGetNumberOfElements(p) != m) )
      throw MatlabException("The initial point for the Lagrange multipliers \
associated with the constraints must be a double-precision array with one \
element for each constraint");
    lambda = mxGetPr(p);
  } else
    lambda = 0;
}
