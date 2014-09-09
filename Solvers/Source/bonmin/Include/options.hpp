// Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         September 25, 2008

// Modified J.Currie September 2011 to suit BONMIN

#ifndef INCLUDE_OPTIONS
#define INCLUDE_OPTIONS

#include "mex.h"
#include "iterate.hpp"
#include "bonminoptions.hpp"
#include "BonTMINLP.hpp"
using namespace  Ipopt;
using namespace Bonmin;

// Class Options.
// -----------------------------------------------------------------
// This class processes the options input from MATLAB.
class Options {
public:

  // The constructor expects as input a point to a MATLAB array, in
  // particular a structure array with the appropriate fields. Note
  // that the Options object does *not* possess an independent copy of
  // some of the MATLAB data (such as the auxiliary data).
  Options (const Iterate& x, Bonmin::BonminSetup& app, 
	   const mxArray* ptr);
  
  // The destructor.
  ~Options();

  // Get the number of variables and the number of constraints.
  friend int numvars        (const Options& options) { return options.n; };
  friend int numconstraints (const Options& options) { return options.m; };
  int numLINconstraints () { return nlin; };
  int numNLconstraints () { return nnlin; };

  // Access the lower and upper bounds on the variables and constraints. 
  const double* lowerbounds () const { return lb; };
  const double* upperbounds () const { return ub; };
  const double* constraintlb() const { return cl; };
  const double* constraintub() const { return cu; };

  // Access the auxiliary data.
  const mxArray* getAuxData() const { return auxdata; };

  // Access the BONMIN options object.
  const BonminOptions bonminOptions() const { return bonmin; };

  // Access the Lagrange multpliers.
  const double* multlb    () const { return zl;     };
  const double* multub    () const { return zu;     };
  const double* multconstr() const { return lambda; };
  
  // Access the MINLP stuff
  const TMINLP::VariableType* variable_types()       const { return var_type; };
  const TNLP::LinearityType*  variable_linearity()   const { return var_lin; };
  const TNLP::LinearityType*  constraint_linearity() const { return cons_lin; };
  
protected:
  int            n;       // The number of optimization variables.
  int            m;       // The number of constraints.
  int            nlin;    // Number of linear constraints
  int            nnlin;   // Number of nonlinear constraints
  double*        lb;      // Lower bounds on the variables.
  double*        ub;      // Upper bounds on the variables.
  double*        cl;      // Lower bounds on constraints.
  double*        cu;      // Upper bounds on constraints.    
  double*        zl;      // Lagrange multipliers for lower bounds.
  double*        zu;      // Lagrange multipliers for upper bounds.
  double*        lambda;  // Lagrange multipliers for constraints.
  const mxArray* auxdata; // MATLAB array containing the auxiliary data.
  BonminOptions  bonmin;   // The BONMIN options.
  
  //MINLP Stuff
  TMINLP::VariableType*  var_type;  // Decision Variable Type (Binary / Integer / Continuous)
  TNLP::LinearityType*   var_lin;   // Decision Variable Linearity
  TNLP::LinearityType*   cons_lin;  // Constraint Variable Linearity

  // These are helper functions used by the class constructor.
  static double* loadLowerBounds      (int n, const mxArray* ptr, 
				       double neginfty);
  static double* loadUpperBounds      (int n, const mxArray* ptr, 
				       double posinfty);
  static int     loadConstraintBounds (const mxArray* ptr, double*& cl, 
				       double*& cu, double neginfty,
				       double posinfty, int &lin, int &nlin);
  static void    loadMultipliers      (int n, int m, const mxArray* ptr, 
				       double*& zl, double*& zu, 
				       double*& lambda);
  
  //MINLP helper functions
  static TMINLP::VariableType* loadVariableTypes(int n, const mxArray *ptr);
  static TNLP::LinearityType*  loadVariableLinearity(int n, const mxArray *ptr);
  static TNLP::LinearityType*  loadConstraintLinearity(int m, int nlin, int nnlin, const mxArray *ptr);
};

#endif
