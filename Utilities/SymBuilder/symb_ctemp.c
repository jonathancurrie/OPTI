/* SYMB_CTEMP - Template for generating SymBuilder C Code Callbacks
 * Copyright (C) 2014 Jonathan Currie (I2C2)                              
 */

#include <mex.h>
#include <ctype.h>
//Prototypes
mwIndex getNoVar();
mwIndex getNoCon();
mwIndex getNNZJac();
mwIndex getNNZHess();
double objective(double *x);
void gradient(double *x, double *v);
void constraints(double *x, double *v);
void jacobian(double *x, double *pr, mwIndex *ir, mwIndex *jc);
void hessian(double *x, double sigma, double *lambda, double *pr, mwIndex *ir, mwIndex *jc);
void printInfo();
void lower(char *str);

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
	double *v, *x, sigma, *lambda, *pr; mwIndex *ir, *jc;
	char *mode; int imode;

    //Check Inputs
    if(nrhs < 1) {
        printInfo();        
        return;
    }   
    if(nrhs < 2) {
        mexErrMsgTxt("You must supply the callback mode and input vector.");
        return;
    }
    if(mxIsEmpty(prhs[0]) || !mxIsChar(prhs[0])) {
        mexErrMsgTxt("The mode must be a string!");
        return;
    }
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
	//Get x
	x = mxGetPr(prhs[1]);
	
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
		plhs[0] = mxCreateSparse(getNoCon(),getNoVar(),getNNZJac(),mxREAL);   
		pr = mxGetPr(plhs[0]);
		ir = mxGetIr(plhs[0]);
		jc = mxGetJc(plhs[0]);
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
		imode = 4;
		plhs[0] = mxCreateSparse(getNoVar(),getNoVar(),getNNZHess(),mxREAL);   
		pr = mxGetPr(plhs[0]);
		ir = mxGetIr(plhs[0]);
		jc = mxGetJc(plhs[0]);
	}
	else
		mexErrMsgTxt("Unknown mode - options are 'obj', 'grad', 'con', 'jac', or 'hess'");
	mxFree(mode);
	
	
	//Call Req Callback
	switch(imode)
	{
		case 0: //objective
			*v = objective(x);
			break;
		case 1: //gradient
			gradient(x,v);
			break;
		case 2: //constraints
			constraints(x,v);
			break;
		case 3: //jacobian
			jacobian(x,pr,ir,jc);
			break;
		case 4: //hessian
			hessian(x,sigma,lambda,pr,ir,jc);
			break;
	}
}

//Print MEX File Info
void printInfo()
{
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" SYMBUILDER C CODE CALLBACK [Built %s]\n",__DATE__);
	mexPrintf("\n Call - symb_ccb(mode,x)\n");
    mexPrintf("\n Modes: 'obj', 'grad', 'con', 'jac' or 'hess'\n");
    mexPrintf("\n J.Currie 2014 (www.i2c2.aut.ac.nz)\n");
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