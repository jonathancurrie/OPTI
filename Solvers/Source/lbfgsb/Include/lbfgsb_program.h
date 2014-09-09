#ifndef INCLUDE_LBFGSB_PROGRAM
#define INCLUDE_LBFGSB_PROGRAM

#include "lbfgsb.h"
#include "program.h"

// Class L-BFGS-B Program.
// -----------------------------------------------------------------
// This is an implementation of the abstract class Program.
class lbfgsb_program: public Program {
public:

    //Constructor
    lbfgsb_program (user_function_data *fun, user_function_data *grad, iter_fun_data *iterF, size_t ndec, 
                    double *lb, double *ub, double *x0, double *x, double *fval, double *iter,
                    int printLevel, int maxIter, int maxFeval, double maxtime, double ftol, int m, double pgtol);

    // The destructor.
    virtual ~lbfgsb_program();

    // These provide definitions for the pure virtual functions of the abstract parent class.
    virtual double computeObjective (int n, double* x);
    virtual void   computeGradient  (int n, double* x, double* g);  
    virtual bool   iterCallback     (int t, double etime, double* x, double f);

    // Run the solver. 
    int runSolver(int &iter, int &feval, double &fval);

protected:
    
    //Inputs
    user_function_data *fun;
    user_function_data *grad;
    iter_fun_data *iterF;
    
    size_t ndec;    //Number of Variables
    
    int noFevals;
    
    int printLevel;
    
    //Outputs
    double *xfinal;
    double *fval;
    double *iter;    
};

#endif
