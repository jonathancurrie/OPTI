// Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         September 25, 2008

// Modified J.Currie September 2011 to suit BONMIN

#ifndef INCLUDE_MATLABPROGRAM
#define INCLUDE_MATLABPROGRAM

#include "iterate.hpp"
#include "options.hpp"
#include "matlabinfo.hpp"
#include "callbackfunctions.hpp"

#include "BonTMINLP.hpp"
using namespace  Ipopt;
using namespace Bonmin;

// Class MatlabProgram
// -----------------------------------------------------------------
class MatlabProgram : public TMINLP {
public:
    
  // The constructor.
  MatlabProgram (const Iterate& x0, const CallbackFunctions& funcs,
		 const Options& options, Iterate& x, const mxArray* auxdata,
		 MatlabInfo& info);
    
  // The destructor.
  virtual ~MatlabProgram();
  
    /** \name Overloaded functions specific to a TMINLP.*/
  //@{
  /** Pass the type of the variables (INTEGER, BINARY, CONTINUOUS) to the optimizer.
     \param n size of var_types (has to be equal to the number of variables in the problem)
  \param var_types types of the variables (has to be filled by function).
  */
  virtual bool get_variables_types(Index n, VariableType* var_types);
 
  /** Pass info about linear and nonlinear variables.*/
  virtual bool get_variables_linearity(Index n, Ipopt::TNLP::LinearityType* var_types);

  /** Pass the type of the constraints (LINEAR, NON_LINEAR) to the optimizer.
  \param m size of const_types (has to be equal to the number of constraints in the problem)
  \param const_types types of the constraints (has to be filled by function).
  */
  virtual bool get_constraints_linearity(Index m, Ipopt::TNLP::LinearityType* const_types);
  
  //NLP Methods
  // Method to return some info about the nonlinear program.
  virtual bool get_nlp_info (int& n, int& m, int& sizeOfJ, int& sizeOfH, 
			     TNLP::IndexStyleEnum& indexStyle);
  
  // Return the bounds for the problem.
  virtual bool get_bounds_info (int n, double* lb, double* ub, int m,
				double* cl, double* cu);
    
  // Return the starting point for the algorithm.
  virtual bool get_starting_point (int n, bool initializeVars, double* vars, 
				   bool initializez, double* zl, double* zu, 
				   int m, bool initializeLambda,
				   double* lambda);
    
  // Compute the value of the objective.
  virtual bool eval_f (int n, const double* vars, bool ignore, double& f);
    
  // Compute the gradient of the objective.
  virtual bool eval_grad_f (int n, const double* vars, bool ignore, 
			    double* grad);
    
  // Evaluate the constraint residuals.
  virtual bool eval_g (int n, const double* vars, bool ignore, int m, 
		       double* g);
    
  // This method either returns: 1.) The structure of the Jacobian
  // (if "Jacobian" is zero), or 2.) The values of the Jacobian (if
  // "Jacobian" is not zero).
  virtual bool eval_jac_g (int numVariables, const double* variables, 
			   bool ignoreThis, int numConstraints, 
			   int sizeOfJ, int* rows, int *cols, double* Jx);
    
  // This method either returns: 1.) the structure of the Hessian of
  // the Lagrangian (if "Hessian" is zero), or 2.) the values of the
  // Hessian of the Lagrangian (if "Hesson" is not zero).
  virtual bool eval_h (int n, const double* vars, bool ignore, double sigma, 
		       int m, const double* lambda, bool ignoretoo, 
		       int sizeOfH, int* rows, int* cols, double* Hx);

  // This method is called when the algorithm is complete.
  virtual void finalize_solution(TMINLP::SolverReturn status,
                                 int n, const double* vars, double obj_value);

  // Intermediate callback method. It is called once per iteration
  // of the IPOPT algorithm.
  virtual bool intermediate_callback (AlgorithmMode mode, int t, double f,
				      double inf_pr, double inf_du,
				      double mu, double d_norm,
				      double regularization_ize,
				      double alpha_du, double alpha_pr,
				      int ls_trials,
				      const IpoptData* ip_data,
				      IpoptCalculatedQuantities* ip_cq);
  
  //Unused
  virtual const SosInfo * sosConstraints() const{return NULL;}
  virtual const BranchingInfo* branchingInfo() const{return NULL;}

protected:
  const Iterate&           x0;       // The initial point.
  const CallbackFunctions& funcs;    // Callback routines.
  const Options&           options;  // Further program info.
  Iterate&                 x;        // Current point.
  const mxArray*           auxdata;  // The auxiliary data.
  MatlabInfo&              info;     // Info passed back to MATLAB.

  // These next two members store information about the structure of
  // the sparse Matlab matrix for the Jacobian of the constraints
  // and the Hessian of the Lagragian.
  SparseMatrix* J;
  SparseMatrix* H;
};

#endif
