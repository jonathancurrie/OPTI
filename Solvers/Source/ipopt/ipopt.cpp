// Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         September 18, 2008

/* Small changes by J.Currie September 2011-2013 */

#include "mex.h"
#include "config_ipopt_default.h"
#include "iterate.hpp"
#include "options.hpp"
#include "matlabexception.hpp"
#include "callbackfunctions.hpp"
#include "matlabinfo.hpp"
#include "matlabjournal.hpp"
#include "matlabprogram.hpp"
#include "IpRegOptions.hpp"
#include "IpJournalist.hpp"
#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"
#include "opti_util.h"
#if defined(LINK_MKL) || defined(LINK_MKL_PARDISO)
    #include "mkl.h"
#endif
#ifdef LINK_MUMPS
    #include "dmumps_c.h"
#endif 

using Ipopt::IsValid;
using Ipopt::RegisteredOption;
using Ipopt::EJournalLevel;
using Ipopt::Journal;
using Ipopt::MatlabJournal;
using Ipopt::IpoptApplication;
using Ipopt::SmartPtr;
using Ipopt::TNLP;
using Ipopt::ApplicationReturnStatus;
using Ipopt::SolveStatistics;

extern void _main();
void printSolverInfo();

// Function definitions.
// -----------------------------------------------------------------
void mexFunction (int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[]) 
  try {

    // Check to see if we have the correct number of input and output arguments.
    if(nrhs < 1) {
        if(nlhs < 1)
            printSolverInfo();
        else
        {
            plhs[0] = mxCreateString(IPOPT_VERSION);  
            plhs[1] = mxCreateDoubleScalar(OPTI_VER);
        }
        return;
    }
    if (nrhs != 3)
      throw MatlabException("Incorrect number of input arguments");
    if (nlhs != 2)
      throw MatlabException("Incorrect number of output arguments");

    // Get the first input which specifies the initial iterate.
    Iterate x0(mxDuplicateArray(prhs[0]));

    // Get the second input which specifies the callback functions.
    CallbackFunctions funcs(prhs[1]);
    
    //Append linear constraints if present
    if(mxGetField(prhs[2],0,"A") && !mxIsEmpty(mxGetField(prhs[2],0,"A"))) 
        funcs.appLinearConstraints(mxGetField(prhs[2],0,"A"));

    // Create a new IPOPT application object and process the options.
    IpoptApplication app(false);
    Options          options(x0,app,prhs[2]);

    // The first output argument is the value of the optimization
    // variables obtained at the solution.
    plhs[0] = mxDuplicateArray(x0);
    Iterate x(plhs[0]);

    // The second output argument stores other information, such as
    // the exit status, the value of the Lagrange multipliers upon
    // termination, the final state of the auxiliary data, and so on.
    MatlabInfo info(plhs[1]);

    // Check to see whether the user provided a callback function for
    // computing the Hessian. This is not needed in the special case
    // when a quasi-Newton approximation to the Hessian is being used.
    if (!options.ipoptOptions().useQuasiNewton() && 
	!funcs.hessianFuncIsAvailable())
      throw MatlabException("You must supply a callback function for computing the Hessian unless you decide to use a quasi-Newton approximation to the Hessian");

    // If the user tried to use her own scaling, report an error.
    if (options.ipoptOptions().userScaling())
      throw MatlabException("The user-defined scaling option does not work in the MATLAB interface for IPOPT");

    // If the user supplied initial values for the Lagrange
    // multipliers, activate the "warm start" option in IPOPT.
    if (options.multlb() && options.multub() &&
	(numconstraints(options)==0 || options.multconstr()) )
      app.Options()->SetStringValue("warm_start_init_point","yes");
     
    //If the user supplied linear constraints, see if we can enable constant jac
    if(funcs.linearAIsAvailable()) {
        int nl = options.numNLconstraints();
        //If no nonlinear constraints, set all linear
        if(!nl) {
            app.Options()->SetStringValue("jac_c_constant","yes");
            app.Options()->SetStringValue("jac_d_constant","yes");
        }
        //Else check if any inequality or equality constraints are nonlinear
        else {
            bool nlEq=false, nlIneq=false;
            const double *cl = options.constraintlb(), *cu = options.constraintub();
            for(int i = 0; i < nl; i++)
                if(cl[i] == cu[i])
                    nlEq = true;
                else
                    nlIneq = true;
            //Enable options based on what we found above
            if(!nlEq) app.Options()->SetStringValue("jac_c_constant","yes");
            if(!nlIneq) app.Options()->SetStringValue("jac_d_constant","yes");
        }
    }
      
    /*** EDIT FOR LIBMWMA57 ***/
    //If using ML MA57, disable MeTiS use (crashes on uncon problem)
    #ifdef LINK_ML_MA57
        if (!options.ipoptOptions().usingMA57()) {
            int value;
            app.Options()->GetIntegerValue("ma57_pivot_order",value,"");
            if(value > 3){
                //mexWarnMsgTxt("overriding ma57 pivot order to option 2 (no MeTiS).");
                app.Options()->SetIntegerValue("ma57_pivot_order",2);}
        }
    #endif
      
    // Set up the IPOPT console.
    EJournalLevel printLevel = (EJournalLevel) 
    options.ipoptOptions().printLevel();
    if(printLevel > 0) { //prevents IPOPT display if we don't want it
        SmartPtr<Journal> console = new MatlabJournal(printLevel);
        app.Jnlst()->AddJournal(console);
    }

    // Intialize the IpoptApplication object and process the options.
    ApplicationReturnStatus exitstatus;
    exitstatus = app.Initialize();
    if (exitstatus != Ipopt::Solve_Succeeded)
      throw MatlabException("IPOPT solver initialization failed");

    // Create a new instance of the constrained, nonlinear program.
    MatlabProgram* matlabProgram 
      = new MatlabProgram(x0,funcs,options,x,options.getAuxData(),info);
    SmartPtr<TNLP> program = matlabProgram;

    // Ask Ipopt to solve the problem.
    exitstatus = app.OptimizeTNLP(program);
    info.setExitStatus(exitstatus);

    // Collect statistics about Ipopt run
    if (IsValid(app.Statistics())) {
      SmartPtr<SolveStatistics> stats = app.Statistics();
      info.setIterationCount(stats->IterationCount());      
      //Get Function Calls
      int obj, con, grad, jac, hess;
      stats->NumberOfEvaluations(obj,con,grad,jac,hess);
      info.setFuncEvals(obj, con, grad, jac, hess);      
      //CPU Time
      info.setCpuTime(stats->TotalCpuTime());
    }

    // Free the dynamically allocated memory.
    mxDestroyArray(x0);

  } catch (std::exception& error) {
    mexErrMsgTxt(error.what());
  }

  
//Print Solver Information
void printSolverInfo()
{    
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" IPOPT: Interior Point Optimizer [v%s]\n",IPOPT_VERSION);
    PRINT_BUILD_INFO;
    mexPrintf("  - Released under the Eclipse Public License: http://opensource.org/licenses/eclipse-1.0\n");
    mexPrintf("  - Source available from: https://projects.coin-or.org/Ipopt\n\n");
    
    mexPrintf(" This binary is statically linked to the following software:\n");
    #ifdef LINK_MUMPS
        mexPrintf("  - MUMPS  [v%s]\n",MUMPS_VERSION);
        mexPrintf("  - METIS  [v4.0.3] (Copyright University of Minnesota)\n");
    #endif
    #ifdef LINK_NETLIB_BLAS
        mexPrintf("  - NETLIB BLAS: http://www.netlib.org/blas/\n  - NETLIB LAPACK: http://www.netlib.org/lapack/\n");
    #endif
    #ifdef LINK_MKL
        mexPrintf("  - Intel Math Kernel Library [v%d.%d R%d]\n",__INTEL_MKL__,__INTEL_MKL_MINOR__,__INTEL_MKL_UPDATE__);        
    #endif
    #ifdef LINK_MKL_PARDISO
        mexPrintf("  - Intel MKL PARDISO [v%d.%d R%d]\n",__INTEL_MKL__,__INTEL_MKL_MINOR__,__INTEL_MKL_UPDATE__);  
    #endif
    #ifdef LINK_MA27
        mexPrintf("  - HSL MA27 (This Binary MUST NOT BE REDISTRIBUTED)\n");
    #endif        
    #ifdef LINK_MA57
        mexPrintf("  - HSL MA57 (This Binary MUST NOT BE REDISTRIBUTED)\n");
        #if defined(LINK_METIS) && !defined(LINK_MUMPS)
            mexPrintf("  - MeTiS [v4.0.3] Copyright University of Minnesota\n");
        #endif
    #endif    
    #if defined(LINK_ML_MA57) || defined(LINK_PARDISO)
        mexPrintf("\n And is dynamically linked to the following software:\n");
        #ifdef LINK_ML_MA57
            mexPrintf("  - MA57    [v3.0] (Included as part of the MATLAB distribution)\n");
        #endif
        #ifdef LINK_PARDISO //BASEL VERSION
            mexPrintf("  - PARDISO [v4.1.2]\n");
        #endif
    #endif
            
    #ifdef HAVE_LINEARSOLVERLOADER
        mexPrintf("\n Dynamically Linked to libhsl.dll containing the HSL Numerical Routines (MA27, MA77, etc)\n   (NOTE: This must exist on the MATLAB path)\n");
    #endif
    
    mexPrintf("\n MEX Interface P.Carbonetto [Modified by J.Currie 2013]\n");
    mexPrintf("-----------------------------------------------------------\n");
}