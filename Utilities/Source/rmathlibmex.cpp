/* RMATHLIBMEX - A MATLAB MEX Interface to the R Math Library
 * Released Under the BSD 3-Clause License:
 * https://www.inverseproblem.co.nz/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2014
 * www.inverseproblem.co.nz
 */

#include <mex.h>
#include <Rmath.h>
#include <string>
#include <ppl.h>
#include "opti_util.h"

using namespace Concurrency;

//Argument Enums (in expected order of arguments)
enum {eFUNC, eIN1, eIN2, eIN3, eIN4};                   
//PRHS Defines    
#define pFUNC    prhs[eFUNC]
#define pIN1     prhs[eIN1]
#define pIN2     prhs[eIN2]
#define pIN3     prhs[eIN3]
#define pIN4     prhs[eIN4]

//Ctrl-C Detection
#ifdef __cplusplus
    extern "C" bool utIsInterruptPending();
    extern "C" void utSetInterruptPending(bool);
#else
    extern bool utIsInterruptPending();
    extern void utSetInterruptPending(bool);
#endif
    
#define max(a,b) (((a) > (b)) ? (a) : (b))

namespace RMathlibMEX
{    
    //Function Prototypes
    void printSolverInfo();
    void checkInputs(const mxArray *prhs[], int nrhs, mwSize &m1, mwSize &n1, mwSize &m2, mwSize &n2, int &fcode);
    void parsolve(double (*fcn)(double), double *in1, double *v, mwSize m1, mwSize n1);
    void parsolve(double (*fcn)(double,double), double *in1, double *in2, double *v, mwSize m1, mwSize n1);
    void parsolve(double (*fcn)(double,double,int), double *in1, double *in2, double *v, mwSize m1, mwSize n1);
    void parsolve(double (*fcn)(double,double,int,int), double *in1, double *in2, double *v, mwSize m1, mwSize n1);
    void parsolve(double (*fcn)(double,double,double,int), double *in1, double *in2, double *in3, double *v, mwSize m1, mwSize n1);
    void parsolve(double (*fcn)(double,double,double,int,int), double *in1, double *in2, double *in3, double *v, mwSize m1, mwSize n1);
    int matchFunc(char *str);
    void lower(char *str);
    
    enum{D1,D2,D2I1,D2I2,D3I1,D3I2};
    //Constants
    int lowerTail = 1;
    int logp = 0;
    
    //Function Structure
    typedef struct {
        char *name;
        int code;
        int args;
        void *pfnc;
        char *detail;
    } RFcn;
    
    RFcn RFcns[] = {
        //T Distribution
        {"\n - T Distribution",-1,-1,NULL," "},
        {"dt",1,D2I1,&dt,"(x,df)   - T Probability Density Function"},
        {"pt",2,D2I2,&pt,"(q,df)   - T Cumulative Distribution Function"},
        {"qt",3,D2I2,&qt,"(p,df)   - T Inverse CDF (Quantile Function)"},
        {"rt",4,D1,&rt,"(df,m,n) - T Distribution Random Number Array"},
        //F Distribution
        {"\n - F Distribution",-1,-1,NULL," "},
        {"df",6,D3I1,&df,"(x,df1,df2)   - F Probability Density Function"},
        {"pf",7,D3I2,&pf,"(q,df1,df2)   - F Cumulative Distribution Function"},
        {"qf",8,D3I2,&qf,"(p,df1,df2)   - F Inverse CDF (Quantile Function)"},
        {"rf",9,D2,&rf,"(df1,df2,m,n) - F Distribution Random Number Array"},
        //Normal Distribution
        {"\n - Normal Distribution",-1,-1,NULL," "},
        {"dnorm",11,D3I1,&dnorm,"(x,mu,sig)   - Normal Probability Density Function"},
        {"pnorm",12,D3I2,&pnorm,"(q,mu,sig)   - Normal Cumulative Distribution Function"},
        {"qnorm",13,D3I2,&qnorm,"(p,mu,sig)   - Normal Inverse CDF (Quantile Function)"},
        {"rnorm",14,D2,&rnorm,"(mu,sig,m,n) - Normal Distribution Random Number Array"},
        //Uniform Distribution
        {"\n - Uniform Distribution",-1,-1,NULL," "},
        {"dunif",16,D3I1,&dunif,"(x,min,max)   - Uniform Probability Density Function"},
        {"punif",17,D3I2,&punif,"(q,min,max)   - Uniform Cumulative Distribution Function"},
        {"qunif",18,D3I2,&qunif,"(p,min,max)   - Uniform Inverse CDF (Quantile Function)"},
        {"runif",19,D2,&runif,"(min,max,m,n) - Uniform Distribution Random Number Array"},
        //Chi^2 Distribution
        {"\n - Chi^2 Distribution",-1,-1,NULL," "},
        {"dchisq",21,D2I1,&dchisq,"(x,df)   - Normal Probability Density Function"},
        {"pchisq",22,D2I2,&pchisq,"(q,df)   - Normal Cumulative Distribution Function"},
        {"qchisq",23,D2I2,&qchisq,"(p,df)   - Normal Inverse CDF (Quantile Function)"},
        {"rchisq",24,D1,&rchisq,"(df,m,n) - Normal Distribution Random Number Array"},
        //Exponential Distribution
        {"\n - Exponential Distribution",-1,-1,NULL," "},
        {"dexp",26,D2I1,&dexp,"(x,rate)   - Exponential Probability Density Function"},
        {"pexp",27,D2I2,&pexp,"(q,rate)   - Exponential Cumulative Distribution Function"},
        {"qexp",28,D2I2,&qexp,"(p,rate)   - Exponential Inverse CDF (Quantile Function)"},
        {"rexp",29,D1,&rexp,"(rate,m,n) - Exponential Distribution Random Number Array"},
        //Beta Distribution
        {"\n - Beta Distribution",-1,-1,NULL," "},
        {"dbeta",31,D3I1,&dbeta,"(x,shape1,shape2)   - Beta Probability Density Function"},
        {"pbeta",32,D3I2,&pbeta,"(q,shape1,shape2)   - Beta Cumulative Distribution Function"},
        {"qbeta",33,D3I2,&qbeta,"(p,shape1,shape2)   - Beta Inverse CDF (Quantile Function)"},
        {"rbeta",34,D2,&rbeta,"(shape1,shape2,m,n) - Beta Distribution Random Number Array"},
        //Cauchy Distribution
        {"\n - Cauchy Distribution",-1,-1,NULL," "},
        {"dcauchy",36,D3I1,&dcauchy,"(x,location,scale)   - Cauchy Probability Density Function"},
        {"pcauchy",37,D3I2,&pcauchy,"(q,location,scale)   - Cauchy Cumulative Distribution Function"},
        {"qcauchy",38,D3I2,&qcauchy,"(p,location,scale)   - Cauchy Inverse CDF (Quantile Function)"},
        {"rcauchy",39,D2,&rcauchy,"(location,scale,m,n) - Cauchy Distribution Random Number Array"},
        //Gamma Distribution
        {"\n - Gamma Distribution",-1,-1,NULL," "},
        {"dgamma",41,D3I1,&dgamma,"(x,shape,scale)   - Gamma Probability Density Function"},
        {"pgamma",42,D3I2,&pgamma,"(q,shape,scale)   - Gamma Cumulative Distribution Function"},
        {"qgamma",43,D3I2,&qgamma,"(p,shape,scale)   - Gamma Inverse CDF (Quantile Function)"},
        {"rgamma",44,D2,&rgamma,"(shape,scale,m,n) - Gamma Distribution Random Number Array"},
        //Log Normal Distribution
        {"\n - Log Normal Distribution",-1,-1,NULL," "},
        {"dlnorm",46,D3I1,&dlnorm,"(x,mu_log,sig_log)   - Log Normal Probability Density Function"},
        {"plnorm",47,D3I2,&plnorm,"(q,mu_log,sig_log)   - Log Normal Cumulative Distribution Function"},
        {"qlnorm",48,D3I2,&qlnorm,"(p,mu_log,sig_log)   - Log Normal Inverse CDF (Quantile Function)"},
        {"rlnorm",49,D2,&rlnorm,"(mu_log,sig_log,m,n) - Log Normal Distribution Random Number Array"},
        //Poisson Distribution
        {"\n - Poisson Distribution",-1,-1,NULL," "},
        {"dpois",51,D2I1,&dpois,"(x,lambda)   - Poisson Probability Density Function"},
        {"ppois",52,D2I2,&ppois,"(q,lambda)   - Poisson Cumulative Distribution Function"},
        {"qpois",53,D2I2,&qpois,"(p,lambda)   - Poisson Inverse CDF (Quantile Function)"},
        {"rpois",54,D1,&rpois,"(lambda,m,n) - Poisson Distribution Random Number Array"},
        //Weibull Distribution
        {"\n - Weibull Distribution",-1,-1,NULL," "},
        {"dweibull",56,D3I1,&dweibull,"(x,shape,scale)   - Weibull Probability Density Function"},
        {"pweibull",57,D3I2,&pweibull,"(q,shape,scale)   - Weibull Cumulative Distribution Function"},
        {"qweibull",58,D3I2,&qweibull,"(p,shape,scale)   - Weibull Inverse CDF (Quantile Function)"},
        {"rweibull",59,D2,&rweibull,"(shape,scale,m,n) - Weibull Distribution Random Number Array"},        
        {NULL,NULL,NULL,NULL,NULL}
    };

    //Main Function
    void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
    {
        //Input Args
        double *p1 = NULL, *p2 = NULL, *p3 = NULL, *t1 = NULL, *t2 = NULL, *t3 = NULL;
        //Output Args
        double *v;
        //Internal Args
        int fcode;
        mwSize i, m1,n1,m2,n2;
        mxArray *p = NULL, *pp = NULL;

        //Check for ver display
        if(nrhs < 1) {
            if(nlhs < 1)
                printSolverInfo();
            else
                plhs[0] = mxCreateString(R_VERSION_STRING);
            return;
        } 

        //Check Inputs
        checkInputs(prhs,nrhs,m1,n1,m2,n2,fcode);        

        //Get Input Args        
        t1 = mxGetPr(pIN1);
        t2 = mxGetPr(pIN2);
        if(nrhs > 3)
            t3 = mxGetPr(pIN3);

        //Check for Scalar Expansion (not random # probs tho)
        if(RFcns[fcode].args != D1 && RFcns[fcode].args != D2) {
            if(m1==1 && n1==1 && (m2 > 1 || n2 > 1)) { //in2 vector/matrix
                p = mxCreateDoubleMatrix(m2,n2,mxREAL); p1 = mxGetPr(p);
                for(i=0;i<m2*n2;i++)
                    p1[i] = *t1;
                p2 = t2;
                p3 = t3;
                m1 = m2; n1 = n2;               
            }
            else if(m2==1 && n2==1 && (m1 > 1 || n1 > 1)) { //in1 vector/matrix
                p = mxCreateDoubleMatrix(m1,n1,mxREAL); p2 = mxGetPr(p);
                for(i=0;i<m1*n1;i++)
                    p2[i] = *t2;
                p1 = t1;
                //May also have to expand p3 here
                if(RFcns[fcode].args == D3I1 || RFcns[fcode].args == D3I2) {
                    pp = mxCreateDoubleMatrix(m1,n1,mxREAL); p3 = mxGetPr(pp);
                    for(i=0;i<m1*n1;i++)
                        p3[i] = *t3;
                }
            }
            else { //no changes required
                p1 = t1;
                p2 = t2;
                p3 = t3;
            }
        }
        else {
            p1 = t1;
            p2 = t2;
            p3 = t3;
        }            
        
        //Create Outputs
        plhs[0] = mxCreateDoubleMatrix(m1,n1, mxREAL);
        v = mxGetPr(plhs[0]);

        //Cast Function Pointer and call corresponding solve routine
        switch(RFcns[fcode].args)
        {
            case D1:   parsolve((double(*)(double))RFcns[fcode].pfnc,p1,v,m1,n1); break;
            case D2:   parsolve((double(*)(double,double))RFcns[fcode].pfnc,p1,p2,v,m1,n1); break;
            case D2I1: parsolve((double(*)(double,double,int))RFcns[fcode].pfnc,p1,p2,v,m1,n1); break;
            case D2I2: parsolve((double(*)(double,double,int,int))RFcns[fcode].pfnc,p1,p2,v,m1,n1); break;
            case D3I1: parsolve((double(*)(double,double,double,int))RFcns[fcode].pfnc,p1,p2,p3,v,m1,n1); break;
            case D3I2: parsolve((double(*)(double,double,double,int,int))RFcns[fcode].pfnc,p1,p2,p3,v,m1,n1); break;
            default:
                mexErrMsgTxt("Not implemented yet");
        }
        
        //Free Memory
        if(p!=NULL) {mxDestroyArray(p); p = NULL;}
        if(pp!=NULL) {mxDestroyArray(pp); pp = NULL;}
    }
    
    //Parallel Solve Routines
    void parsolve(double (*fcn)(double), double *in1, double *v, mwSize m1, mwSize n1)
    {
        //Scalar
        if(m1 == 1 && n1 == 1)
        {
            v[0] = (*fcn)(in1[0]);
            return;
        }        
        //Vector
        if((m1 > 1 && n1 == 1) || (n1 > 1 && m1 == 1))
        {
            size_t m = max(m1,n1);
            parallel_for((size_t)0,m,(size_t)1,[&](size_t i)
            {     
                v[i] = (*fcn)(in1[0]);
            });
            return;
        }
        //Matrix
        parallel_for((size_t)0,m1,(size_t)1,[&](size_t i)
        {
            for(size_t j = 0; j < n1; j++)
            {
                size_t ind = i+j*m1;
                v[ind] = (*fcn)(in1[0]);
            }   
        });
    } 
    void parsolve(double (*fcn)(double,double), double *in1, double *in2, double *v, mwSize m1, mwSize n1)
    {
        //Scalar
        if(m1 == 1 && n1 == 1)
        {
            v[0] = (*fcn)(in1[0],in2[0]);
            return;
        }        
        //Vector
        if((m1 > 1 && n1 == 1) || (n1 > 1 && m1 == 1))
        {
            size_t m = max(m1,n1);
            parallel_for((size_t)0,m,(size_t)1,[&](size_t i)
            {     
                v[i] = (*fcn)(in1[0],in2[0]);
            });
            return;
        }
        //Matrix
        parallel_for((size_t)0,m1,(size_t)1,[&](size_t i)
        {
            for(size_t j = 0; j < n1; j++)
            {
                size_t ind = i+j*m1;
                v[ind] = (*fcn)(in1[0],in2[0]);
            }   
        });
    } 
    
    void parsolve(double (*fcn)(double,double,int), double *in1, double *in2, double *v, mwSize m1, mwSize n1)
    {
        //Scalar
        if(m1 == 1 && n1 == 1)
        {
            v[0] = (*fcn)(in1[0],in2[0],logp);
            return;
        }        
        //Vector
        if((m1 > 1 && n1 == 1) || (n1 > 1 && m1 == 1))
        {
            size_t m = max(m1,n1);
            parallel_for((size_t)0,m,(size_t)1,[&](size_t i)
            {     
                v[i] = (*fcn)(in1[i],in2[i],logp);
            });
            return;
        }
        //Matrix
        parallel_for((size_t)0,m1,(size_t)1,[&](size_t i)
        {
            for(size_t j = 0; j < n1; j++)
            {
                size_t ind = i+j*m1;
                v[ind] = (*fcn)(in1[ind],in2[ind],logp);
            }   
        });
    } 
    
    void parsolve(double (*fcn)(double,double,int,int), double *in1, double *in2, double *v, mwSize m1, mwSize n1)
    {
        //Scalar
        if(m1 == 1 && n1 == 1)
        {
            v[0] = (*fcn)(in1[0],in2[0],lowerTail,logp);
            return;
        }        
        //Vector
        if((m1 > 1 && n1 == 1) || (n1 > 1 && m1 == 1))
        {
            size_t m = max(m1,n1);
            parallel_for((size_t)0,m,(size_t)1,[&](size_t i)
            {     
                v[i] = (*fcn)(in1[i],in2[i],lowerTail,logp);
            });
            return;
        }
        //Matrix
        parallel_for((size_t)0,m1,(size_t)1,[&](size_t i)
        {
            for(size_t j = 0; j < n1; j++)
            {
                size_t ind = i+j*m1;
                v[ind] = (*fcn)(in1[ind],in2[ind],lowerTail,logp);
            }   
        });
    }   
    
    void parsolve(double (*fcn)(double,double,double,int), double *in1, double *in2, double *in3, double *v, mwSize m1, mwSize n1)
    {
        //Scalar
        if(m1 == 1 && n1 == 1)
        {
            v[0] = (*fcn)(in1[0],in2[0],in3[0],logp);
            return;
        }        
        //Vector
        if((m1 > 1 && n1 == 1) || (n1 > 1 && m1 == 1))
        {
            size_t m = max(m1,n1);
            parallel_for((size_t)0,m,(size_t)1,[&](size_t i)
            {     
                v[i] = (*fcn)(in1[i],in2[i],in3[i],logp);
            });
            return;
        }
        //Matrix
        parallel_for((size_t)0,m1,(size_t)1,[&](size_t i)
        {
            for(size_t j = 0; j < n1; j++)
            {
                size_t ind = i+j*m1;
                v[ind] = (*fcn)(in1[ind],in2[ind],in3[ind],logp);
            }   
        });
    } 
    
    void parsolve(double (*fcn)(double,double,double,int,int), double *in1, double *in2, double *in3, double *v, mwSize m1, mwSize n1)
    {
        //Scalar
        if(m1 == 1 && n1 == 1)
        {
            v[0] = (*fcn)(in1[0],in2[0],in3[0],lowerTail,logp);
            return;
        }        
        //Vector
        if((m1 > 1 && n1 == 1) || (n1 > 1 && m1 == 1))
        {
            size_t m = max(m1,n1);
            parallel_for((size_t)0,m,(size_t)1,[&](size_t i)
            {     
                v[i] = (*fcn)(in1[i],in2[i],in3[i],lowerTail,logp);
            });
            return;
        }
        //Matrix
        parallel_for((size_t)0,m1,(size_t)1,[&](size_t i)
        {
            for(size_t j = 0; j < n1; j++)
            {
                size_t ind = i+j*m1;
                v[ind] = (*fcn)(in1[ind],in2[ind],in3[ind],lowerTail,logp);
            }   
        });
    } 

    void checkInputs(const mxArray *prhs[], int nrhs, mwSize &m1, mwSize &n1, mwSize &m2, mwSize &n2, int &fcode)
    {   
        char *func;
        if(nrhs < 3)
            mexErrMsgTxt("You must supply at least 3 arguments to rmathlib!\n\rmathlib('fcn',in1,in2)\n");                
        //Check data type
        if(mxIsEmpty(pFUNC) || !mxIsChar(pFUNC))
            mexErrMsgTxt("Function name must be a string!");
        if(mxIsEmpty(pIN1) || !mxIsDouble(pIN1) || mxIsSparse(pIN1) || mxIsComplex(pIN1))
            mexErrMsgTxt("Input 1 must be a real, dense, double scalar/vector/maxtrix");
        if(mxIsEmpty(pIN2) || !mxIsDouble(pIN2) || mxIsSparse(pIN2) || mxIsComplex(pIN2))
            mexErrMsgTxt("Input 2 must be a real, dense, double scalar/vector/maxtrix");                
        
        //Get Sizes
        m1 = mxGetM(pIN1); n1 = mxGetN(pIN1);
        m2 = mxGetM(pIN2); n2 = mxGetN(pIN2);
        
        //Find and match our function
        func = mxArrayToString(pFUNC);
        fcode = matchFunc(func);         
        mxFree(func);
        
        //Depending on function found, check for correct number of args and sizes
        if(RFcns[fcode].args == D1) { //random no gen e.g. rt
            if(nrhs < 4)
                mexErrMsgTxt("This function requires at least four inputs - rmathlib('fun',df,m,n)");
            if(mxGetNumberOfElements(pIN1) > 1)
                mexErrMsgTxt("Input 1 must be a real, dense, double scalar (degrees of freedom)");
            if(mxGetNumberOfElements(pIN2) > 1)
                mexErrMsgTxt("Input 2 must be a real, dense, double scalar (# rows)");
            if(mxIsEmpty(pIN3) || !mxIsDouble(pIN3) || mxIsSparse(pIN3) || mxIsComplex(pIN3) || mxGetNumberOfElements(pIN3) > 1)
                mexErrMsgTxt("Input 3 must be a real, dense, double scalar (# cols)"); 
            //Update Sizes
           m1 = (mwSize)*mxGetPr(pIN2);
           n1 = (mwSize)*mxGetPr(pIN3);
        }
        else if(RFcns[fcode].args == D2) { //random no gen e.g. rf
            if(nrhs < 5)
                mexErrMsgTxt("This function requires at least five inputs - (e.g.) rmathlib('fun',df1,df2,m,n)");
            if(mxGetNumberOfElements(pIN1) > 1)
                mexErrMsgTxt("Input 1 must be a real, dense, double scalar (degrees of freedom 1)");
            if(mxGetNumberOfElements(pIN1) > 2)
                mexErrMsgTxt("Input 2 must be a real, dense, double scalar (degrees of freedom 2)");
            if(mxIsEmpty(pIN3) || !mxIsDouble(pIN3) || mxIsSparse(pIN3) || mxIsComplex(pIN3) || mxGetNumberOfElements(pIN3) > 1)
                mexErrMsgTxt("Input 3 must be a real, dense, double scalar (# rows)");
            if(mxIsEmpty(pIN4) || !mxIsDouble(pIN4) || mxIsSparse(pIN4) || mxIsComplex(pIN4) || mxGetNumberOfElements(pIN4) > 1)
                mexErrMsgTxt("Input 4 must be a real, dense, double scalar (# cols)"); 
            //Update Sizes
           m1 = (mwSize)*mxGetPr(pIN3);
           n1 = (mwSize)*mxGetPr(pIN4);
        }
        else if(RFcns[fcode].args == D2I1 || RFcns[fcode].args == D2I2 || RFcns[fcode].args == D3I1 || RFcns[fcode].args == D3I2) {            
            //OK if both scalar, one scalar, or both vector/matrix of same size
            if((m1 > 1 || n1 > 1) && (m2 > 1 || n2 > 1)) {
                if(m1*n1 != m2*n2)
                    mexErrMsgTxt("If both inputs are non-scalar, they must have the same number of elements");
                if(m1!=m2 || n1!=n2)
                    mexErrMsgTxt("If both inputs are non-scalar, they must have the same size dimensions");
            }         
            //Extra check for D3 functions
            if(RFcns[fcode].args == D3I1 || RFcns[fcode].args == D3I2) {
                if(nrhs < 4)
                    mexErrMsgTxt("This function requires at least four inputs - (e.g.) rmathlib('fun',x,df1,df2)");
                if(mxIsEmpty(pIN3) || !mxIsDouble(pIN3) || mxIsSparse(pIN3) || mxIsComplex(pIN3))
                    mexErrMsgTxt("Input 3 must be a real, dense, double scalar/vector/maxtrix");                 
                //Get New Sizes
                mwSize m3 = mxGetM(pIN3), n3 = mxGetN(pIN3);
                //DF1 must have the same dimensions as DF2
                if(m2 != m3 || n2 != n3)
                    mexErrMsgTxt("Both inputs 2 and 3 (e.g.) df1,df2 must have the same dimensions");
            }
        }
        else {
            mexPrintf("Func Code: %d, Args: %d\n",fcode,RFcns[fcode].args);
            mexErrMsgTxt("Unknown function form (internal error)");
        }
    }

    //Print Solver Information
    void printSolverInfo()
    {    
        char vbuf[6]; getVSVer(vbuf);  
        mexPrintf("\n-----------------------------------------------------------\n");
        mexPrintf(" RMATHLIB: R Math Library [v%s, Built %s, VS%s]\n",R_VERSION_STRING,__DATE__,vbuf);              
        mexPrintf("  - Released under the GNU General Public License: http://www.gnu.org/licenses/gpl-2.0.html\n");
        mexPrintf("  - Source available from: http://cran.stat.sfu.ca/\n");

        mexPrintf("\n Example Use: rmathlib('qt',0.95,1:10)\n");
        
        mexPrintf("\n Available Functions:\n");
        RFcn *r = &RFcns[0];
        //Match the name
        while(r && r->name) {
            mexPrintf("%s%s\n",r->name,r->detail);
            r++;
        }
        
        mexPrintf("\n where\n  df  = Degrees of Freedom\n");
        mexPrintf("  m   = Number of Rows\n");
        mexPrintf("  n   = Number of Columns\n");
        mexPrintf("  mu  = Mean\n");
        mexPrintf("  sig = Standard Deviation\n");
        mexPrintf("\n Note all functions are vectorized and parallelized, and return the lower tail by default.\n");
        
        mexPrintf("\n MEX Interface J.Currie 2014 [BSD3] (www.inverseproblem.co.nz)\n");
        mexPrintf("-----------------------------------------------------------\n");
    }
    
    //Match Input Function Name
    int matchFunc(char *str)
    {
        //Get First Entry
        RFcn *r = &RFcns[0];
        //Lowercase entry
        lower(str);
        //Match the name
        while(r && r->name)
        {
            if(strcmp(str, r->name) == 0)
                return r->code;
            r++;
        }         
        mexErrMsgTxt("Unknown RMathlib function! Use rmathlib() to see all available functions.");
    }
    
    //Convert input string to lowercase
    void lower(char *str)
    {
        int i = 0;
        while(str[i]) {
            str[i] = tolower(str[i]);  
            i++;
        }
    }
}
