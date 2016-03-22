/* SYMB_CTEMP - Template for generating SymBuilder C Code Callbacks
 * This Version also Uses CppAD for generating Derivatives
 * Copyright (C) 2014-2016 Jonathan Currie (I2C2)                              
 */

#include <mex.h>
#include <string.h>
#include <math.h>
#include <cppad/cppad.hpp>

using std::vector;
using std::set;

//Prototypes
mwIndex getNoVar();
mwIndex getNoCon();
template <typename Type> Type objective(const vector<Type> &x);
template <typename Type> void constraints(const vector<Type> &x, vector<Type> &v);
void printInfo();
void lower(char *str);

//Global Vars (to keep memory between calls)
static vector<double> xvec, cvec, jac, hes, w;
static vector<size_t> jrow, jcol, hrow, hcol;
static vector<set<size_t>> jstr, hstr;
static CppAD::ADFun<double> obj,con,lag;
static CppAD::sparse_jacobian_work jwork;
static CppAD::sparse_hessian_work hwork;
static size_t nnzJac=0, nnzHess=0, nnzHessLT=0;
static mwIndex *jir=NULL, *jjc=NULL, *hir=NULL, *hjc=NULL;

//Memory Clear Function
static void mexExit(void)
{
    //mexPrintf("Clearing memory...\n");
    if(!xvec.empty())
        xvec.~vector();
    if(!cvec.empty())
        cvec.~vector();
    if(jir!=NULL) {
        mxFree(jir);
        jir=NULL;
    }
    if(jjc!=NULL) {
        mxFree(jjc);
        jjc=NULL;
    }
    if(!jrow.empty())
        jrow.~vector();
    if(!jcol.empty())
        jcol.~vector();
    if(!jac.empty())
        jac.~vector();
    if(hir!=NULL) {
        mxFree(hir);
        hir=NULL;
    }
    if(hjc!=NULL) {
        mxFree(hjc);
        hjc=NULL;
    }
    if(!hrow.empty())
        hrow.~vector();
    if(!hcol.empty())
        hcol.~vector();
    if(!hes.empty())
        hes.~vector();
    if(!w.empty())
        w.~vector();
}

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
	double *v, *x, sigma, *lambda, *pr;
	char *mode; int imode;

    //Check Inputs
    if(nrhs < 1) {
        printInfo();        
        return;
    }   
    
    if(mxIsEmpty(prhs[0]) || !mxIsChar(prhs[0])) {
        mexErrMsgTxt("The mode must be a string!");
        return;
    }
    
    //If we have x, check it
    if(nrhs > 1) {
        if(!mxIsEmpty(prhs[1])) {
            if(mxIsClass(prhs[1],"scipvar") || mxIsClass(prhs[1],"barvec")) {
                mexErrMsgTxt("SCIP and BARON cannot be used with this callback function - please specify 'mcode' via symbset as the cbmode.");
                return;
            }           
            if(!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1])) {
                mexErrMsgTxt("The input vector must be a dense real double vector!");
                return;
            }
        }
        else {
            mexErrMsgTxt("The input vector must be a dense real double vector!");
            return;
        }
        //Check x input size
        if(mxGetNumberOfElements(prhs[1]) != getNoVar()) {
            mexErrMsgTxt("The input vector is not the right size!");
        }
        //Allocate memory, if required
        if(xvec.empty())
            xvec.resize(getNoVar());        
        //Get x and copy to xvec
        x = mxGetPr(prhs[1]);
        memcpy(&xvec[0],x,getNoVar()*sizeof(double));
    }
	
	//Determine input mode and setup return variable
	mode = mxArrayToString(prhs[0]);
	lower(mode);
	if(!strcmp(mode,"obj")) {
		imode = 0;
		plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
		v = mxGetPr(plhs[0]);
	}
	else if(!strcmp(mode,"grad")) {
		imode = 1;
		plhs[0] = mxCreateDoubleMatrix(1,getNoVar(), mxREAL);
		v = mxGetPr(plhs[0]);
	}
	else if(!strcmp(mode,"con")) {
		imode = 2;
		plhs[0] = mxCreateDoubleMatrix(getNoCon(),1, mxREAL);
		v = mxGetPr(plhs[0]);
	}
	else if(!strcmp(mode,"jac")) {
		imode = 3;
        //Can't allocate here until we know sparsity pattern		
	}
    else if(!strcmp(mode,"jacstr")) {
		imode = 4;
        //Can't allocate here until we know sparsity pattern		
	}
    else if(!strcmp(mode,"hess")) {
		if(nrhs < 4) {
			mexErrMsgTxt("You must supply the callback mode, input vector, sigma and lambda for Hessian Evaluations.");
			return;
		}
		//Check length of Sigma
		if(mxIsEmpty(prhs[2]) || !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 1)
			mexErrMsgTxt("Sigma must be a real, double scalar.");
		//Check length of Lambda
		if(!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxIsSparse(prhs[3]) || mxGetNumberOfElements(prhs[3]) != getNoCon())
			mexErrMsgTxt("Lambda must be a real, double, dense vector with ncon elements.");
		//Get Sigma, Lambda
		sigma = *mxGetPr(prhs[2]);
        lambda = mxGetPr(prhs[3]);
		imode = 5;		
        //Can't allocate here until we know sparsity pattern	
	}
    else if(!strcmp(mode,"hstr")) {
        imode = 6;
        //Can't allocate here until we know sparsity pattern	
    }
	else
		mexErrMsgTxt("Unknown mode - options are 'obj', 'grad', 'con', 'jac', 'jacstr', 'hess' or 'hstr'");
	mxFree(mode);
	
	//Ensure we did have x for normal callbacks
    if(imode != 4 && imode != 6 && nrhs < 2)
        mexErrMsgTxt("You must supply the callback mode and input vector.");
    
	//Call Req Callback
	switch(imode)
	{
		case 0: //objective            
			*v = objective(xvec);
			break;
		case 1: //gradient
            //Check if we have recorded the objective yet
            if(obj.Memory()==0) { //new, tape operations
                vector< CppAD::AD<double> > X(getNoVar());
                memcpy(&X[0],x,getNoVar()*sizeof(double));
                CppAD::Independent(X);
                vector< CppAD::AD<double> > Y(1);
                Y[0] = objective(X);     
                obj = CppAD::ADFun<double>(X, Y);
                //obj.optimize();
                mexAtExit(mexExit); //also register memory clear function
                //mexPrintf("Evaluated Tape for Gradient\n");
            }
            //Evaluate "Jacobian" for gradient
            memcpy(v,&(obj.Jacobian(xvec)[0]),getNoVar()*sizeof(double));
			break;
		case 2: //constraints
            //Check if we have constraint memory yet
            if(cvec.empty())
                cvec.resize(getNoCon()); //allocate it
            //Evaluate Constraints
			constraints(xvec,cvec);
            //Copy Out
            memcpy(v,&cvec[0],getNoCon()*sizeof(double));
			break;
		case 3: //jacobian
        case 4: //jacobian structure
			//Check if we have recorded the constraints yet
            if(con.Memory()==0){ //new, tape operations
                vector< CppAD::AD<double> > X(getNoVar());
                memcpy(&X[0],x,getNoVar()*sizeof(double));
                CppAD::Independent(X);
                vector< CppAD::AD<double> > Y(getNoCon());
                constraints(X,Y);     
                con = CppAD::ADFun<double>(X, Y);
                //con.optimize();
                mexAtExit(mexExit); //also register memory clear function
                //mexPrintf("Evaluated Tape for Jacobian\n");
            }
            //Check if we have the sparsity pattern yet
            if(jstr.empty()) {                
                vector<set<size_t>> r(getNoVar());
                for(size_t i = 0; i < getNoVar(); i++)
                    r[i].insert(i); //identity matrix 
                jstr.resize(getNoCon());
                jstr = con.ForSparseJac(getNoVar(),r,true); //note transpose
                //Determine nnzs
                for(int i = 0; i < jstr.size(); i++)
                    nnzJac += jstr[i].size();
                //Save ir, jc for jac
                jir = (mwIndex*)mxCalloc(nnzJac,sizeof(mwIndex));
                jjc = (mwIndex*)mxCalloc(getNoVar()+1,sizeof(mwIndex));                
                mexMakeMemoryPersistent(jir);
                mexMakeMemoryPersistent(jjc);
                jwork.clear(); //reset jacobian calculations
                //Col starts
                jjc[0] = 0;
                for(int i = 1; i <= getNoVar(); i++)
                    jjc[i] = (mwIndex)(jjc[i-1] + jstr[i-1].size());
                //Rows
                size_t idx = 0;
                for(int i = 0; i < jstr.size(); i++)
                    for (set<size_t>::iterator it=jstr[i].begin(); it!=jstr[i].end(); ++it)
                        jir[idx++] = (mwIndex)*it;
                //Build missing triple so we can eval just sparse elements of Jac
                jrow.resize(nnzJac);
                jcol.resize(nnzJac);
                idx = 0;
                for(size_t i = 0; i < nnzJac; i++)
                    jrow[i] = jir[i];
                for(size_t i = 0; i < getNoVar();i++)
                    for(size_t j = jjc[i]; j < jjc[i+1]; j++)
                        jcol[idx++] = i;
                //Re-do with no transpose... (bad really...)
                jstr = con.ForSparseJac(getNoVar(),r,false); 
                //mexPrintf("Determined Jac Sparsity Structure (%d nzs)\n",nnzJac);
            }
            //Create Sparse Return Matrix
            plhs[0] = mxCreateSparse(getNoCon(),getNoVar(),nnzJac,mxREAL);   
            pr = mxGetPr(plhs[0]);
            memcpy(mxGetIr(plhs[0]),jir,nnzJac*sizeof(mwIndex));
            memcpy(mxGetJc(plhs[0]),jjc,(getNoVar()+1)*sizeof(mwIndex));
            //If we want the sparsity pattern only, fill in return matrix with 1s
            if(imode==4) {   
                for(int i = 0; i < nnzJac; i++)
                	pr[i] = 1.0;                
            }
            //Else, evaluate sparse jacobian and return as sparse matrix
            else {
                //Check if we have jacobian memory yet
                if(jac.empty())
                    jac.resize(nnzJac); //allocate it
                //If ndec > ncon, use reverse mode
                if(getNoVar() > getNoCon())                   
                    con.SparseJacobianReverse(xvec,jstr,jrow,jcol,jac,jwork);
                //else use forward
                else
                    con.SparseJacobianForward(xvec,jstr,jrow,jcol,jac,jwork);
                //Copy out
                memcpy(pr,&jac[0],nnzJac*sizeof(double));
            }
			break;
		case 5: //hessian of the lagrangian
        case 6: //hessian structure
            //Check if we have recorded the objective+constraints yet
            //Not sure if we can reuse ones we have done above??
            if(lag.Memory()==0){ //new, tape operations
                vector< CppAD::AD<double> > X(getNoVar());
                memcpy(&X[0],x,getNoVar()*sizeof(double));
                CppAD::Independent(X);
                //Output Array
                vector< CppAD::AD<double> > Y(1); 
                vector< CppAD::AD<double> > Yc(getNoCon()); 
                Y[0] = objective(X); //eval objective   
                if(getNoCon() > 0)
                    constraints(X,Yc); //eval constraints     
                Yc.insert(Yc.begin(),Y.begin(),Y.end());
                //Create ADFun
                lag.Dependent(Yc);
                //lag.optimize();
                mexAtExit(mexExit); //also register memory clear function
                //mexPrintf("Evaluated Tape for Hessian\n");
            }
            //Check if we have the sparsity pattern yet
            if(hstr.empty()) {        
                //First eval jac structure (not sure why)
                vector< std::set<size_t> > r(getNoVar());
                for(size_t i = 0; i < getNoVar(); i++)
                    r[i].insert(i);
                lag.ForSparseJac(getNoVar(), r);
                //Now do Hessian structure
                vector<set<size_t>> s(1);
                for(size_t i = 0; i < getNoCon()+1; i++)
                    s[0].insert(i); //identity matrix 
                hstr.resize(getNoVar());
                hstr = lag.RevSparseHes(getNoVar(),s); 
                //Determine total nnzs
                for(int i = 0; i < hstr.size(); i++)
                    nnzHess += hstr[i].size();
                //Determine nnzs in lower tri
                for(int i = 0; i < hstr.size(); i++)
                    for (set<size_t>::iterator it=hstr[i].begin(); it!=hstr[i].end(); ++it)
                        if(*it >= i)
                            nnzHessLT++;
                
                //Save ir, jc for jac
                hir = (mwIndex*)mxCalloc(nnzHessLT,sizeof(mwIndex));
                hjc = (mwIndex*)mxCalloc(getNoVar()+1,sizeof(mwIndex));                
                mexMakeMemoryPersistent(hir);
                mexMakeMemoryPersistent(hjc);
                hwork.clear(); //reset hessian calculations
                //Col & Row Starts
                size_t idx = 0;
                for(int i = 0; i < hstr.size(); i++) {
                    hjc[i] = idx;
                    for (set<size_t>::iterator it=hstr[i].begin(); it!=hstr[i].end(); ++it)
                        if(*it >= i)
                            hir[idx++] = (mwIndex)*it;
                }
                hjc[getNoVar()] = nnzHessLT;
                //Build missing triple so we can eval just sparse elements of Jac
                hrow.resize(nnzHessLT);
                hcol.resize(nnzHessLT);
                idx = 0;
                for(size_t i = 0; i < nnzHessLT; i++)
                    hrow[i] = hir[i];
                for(size_t i = 0; i < getNoVar();i++)
                    for(size_t j = hjc[i]; j < hjc[i+1]; j++)
                        hcol[idx++] = i;
                //mexPrintf("Determined Hess Sparsity Structure (%d nzs in tril)\n",nnzHessLT);
            }
            //Create Sparse Return Matrix
            plhs[0] = mxCreateSparse(getNoVar(),getNoVar(),nnzHessLT,mxREAL);   
            pr = mxGetPr(plhs[0]);
            memcpy(mxGetIr(plhs[0]),hir,nnzHessLT*sizeof(mwIndex));
            memcpy(mxGetJc(plhs[0]),hjc,(getNoVar()+1)*sizeof(mwIndex));
            //If we want the sparsity pattern only, fill in return matrix with 1s
            if(imode==6) {   
                for(int i = 0; i < nnzHessLT; i++)
                	pr[i] = 1.0;                
            }
            //Else, evaluate sparse hessian and return as sparse matrix
            else {
                //Check if we have hessian memory yet
                if(hes.empty())
                    hes.resize(nnzHessLT); //allocate it  
                if(w.empty())
                    w.resize(1+getNoCon()); //allocate it
                //Copy in Weights
                w[0] = sigma;
                for(int i = 0; i < getNoCon(); i++)
                    w[i+1] = lambda[i];
                //If ndec > ncon, use reverse mode
                lag.SparseHessian(xvec,w,hstr,hrow,hcol,hes,hwork);
                //Copy out elements
                memcpy(pr,&hes[0],nnzHessLT*sizeof(double));
            }
			break;
	}
}

//Print MEX File Info
void printInfo()
{
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" SYMBUILDER C CODE CALLBACK with CppAD [%s] [Built %s]\n",CPPAD_PACKAGE_STRING,__DATE__);
    mexPrintf("  - CppAD Released under the Eclipse Public License: http://opensource.org/licenses/eclipse-1.0\n");
    mexPrintf("  - CpPAD Source available from: http://www.coin-or.org/CppAD\n");
	mexPrintf("\n Call - symb_ccb(mode,x)\n");
    mexPrintf("\n Modes: 'obj', 'grad', 'con', 'jac' 'jacstr', 'hess' or 'hstr'\n");
    mexPrintf("\n J.Currie 2014-2016 (www.i2c2.aut.ac.nz)\n");
    mexPrintf("-----------------------------------------------------------\n");
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