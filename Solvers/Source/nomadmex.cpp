/*-------------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct search MATLAB Interface     */
/*                                                                                     */
/*  Copyright (C) 2001-2013  Mark Abramson        - the Boeing Company, Seattle        */
/*                           Charles Audet        - Ecole Polytechnique, Montreal      */
/*                           Gilles Couture       - Ecole Polytechnique, Montreal      */
/*                           John Dennis          - Rice University, Houston           */
/*                           Sebastien Le Digabel - Ecole Polytechnique, Montreal      */
/*                           Christophe Tribes    - Ecole Polytechnique, Montreal      */
/*                                                                                     */
/*  funded in part by AFOSR and Exxon Mobil                                            */
/*                                                                                     */
/*  MEX Interface Author: Jonathan Currie 2012-2013 (I2C2)                             */
/*                                                                                     */
/*  MEX Interface Contact information:                                                 */
/*    jocurrie@aut.ac.nz                                                               */
/*                                                                                     */
/*  This program is free software: you can redistribute it and/or modify it under the  */
/*  terms of the GNU Lesser General Public License as published by the Free Software   */
/*  Foundation, either version 3 of the License, or (at your option) any later         */
/*  version.                                                                           */
/*                                                                                     */
/*  This program is distributed in the hope that it will be useful, but WITHOUT ANY    */
/*  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A    */
/*  PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.   */
/*                                                                                     */
/*  You should have received a copy of the GNU Lesser General Public License along     */
/*  with this program. If not, see <http://www.gnu.org/licenses/>.                     */
/*                                                                                     */
/*  You can find information on the NOMAD software at www.gerad.ca/nomad               */
/*-------------------------------------------------------------------------------------*/

#define NOMADMEX_VERSION "1.23  [23/06/14]"
//NOTE from Version 1.15 on this MEX file contains a dual interface:

// 1) The default (GERAD) interface is:
//      [x,fval,exitflag,iter] = nomad(fun,x0,lb,ub,opts)

// 2) Otherwise if OPTI_VERSION is defined it will compile the OPTI Toolbox interface:
//      [x,fval,exitflag,iter] = nomad(fun,x0,lb,ub,nlcon,nlrhs,xtype,opts)

// Both contain near identical functionality, except the GERAD version includes nonlinear
// constraints with the objective (blackbox) function. The OPTI version separates them,
// and allows a user specified rhs. The OPTI version also rounds x0 for integer constraints.

#include "mex.h"
#include "nomad.hpp"
#include "defines.hpp"
#include <stdio.h>
#include <string.h>

using namespace std;

//Function handle structure
#define FLEN 128 /* max length of user function name */
#define MAXRHS 3 /* max nrhs for user function */
typedef struct {
     char f[FLEN];
     mxArray *plhs[1];
     mxArray *prhs[MAXRHS];
     int xrhs, nrhs;
     double *nlrhs;
} usrFcn;

//Iteration callback structure
typedef struct {
    char f[FLEN];
    mxArray *plhs[1];
    mxArray *prhs[4];
    bool enabled;
} iter_fun_data;

//Ctrl-C Detection
#ifdef __cplusplus
    extern "C" bool utIsInterruptPending();
    extern "C" void utSetInterruptPending(bool);
#else
    extern bool utIsInterruptPending();
    extern void utSetInterruptPending(bool);
#endif
    
//Argument Enums (in expected order of arguments)
#ifdef OPTI_VERSION
    #include "opti_util.h"
    enum {eFUN, eX0, eLB, eUB, eNLCON, eNLRHS, eXTYPE, eOPTS};
#else //GERAD VERSION
    enum {eFUN, eX0, eLB, eUB, eOPTS};
    enum {eNLCON, eNLRHS, eXTYPE}; //placeholders
#endif
//PRHS Defines    
#define pFUN    prhs[eFUN]
#define pX0     prhs[eX0]
#define pLB     prhs[eLB]
#define pUB     prhs[eUB]
#define pNLCON  prhs[eNLCON]
#define pNLRHS  prhs[eNLRHS]
#define pXTYPE  prhs[eXTYPE]
#define pOPTS   prhs[eOPTS]

//Function Prototypes
void printSolverInfo();
int checkInputs(const mxArray *prhs[], int nrhs, mxArray *plhs[], int nlhs);
vector<NOMAD::bb_input_type> detInTypes(char *xtype, size_t n);
vector<NOMAD::bb_output_type> detRetTypes(usrFcn *fun, mxArray *out_types, int *nobj, usrFcn *con, int *ncon, double *x0, size_t n);
void setNOMADopts(NOMAD::Parameters &p, const mxArray *opts);
bool isNMDOption(const char *field);
NOMAD::bb_output_type getOutEnum(char *str);
void lower(char *str);
double getStatus(int stat);

//MATLAB Evaluator Class
class matlabEval : public NOMAD::Evaluator {
private:
    usrFcn *fun;
    bool hasSur;
    int nobj;
    usrFcn *con;
    int ncon;
    iter_fun_data *iterF;
    int citer;
    
public:
    //Constructor
    matlabEval(const NOMAD::Parameters &p, usrFcn *_obj, int _nobj, usrFcn *_con, int _ncon, iter_fun_data *_iterF) : NOMAD::Evaluator(p)
    {
        fun     = _obj;
        hasSur  = p.has_sgte();
        nobj    = _nobj;
        con     = _con;   
        ncon    = _ncon;
        iterF   = _iterF;
        citer   = 0;
    }        
    //Destructor
    ~matlabEval(void) {}
    //Function + Constraint Evaluation
    bool eval_x(NOMAD::Eval_Point &x, const NOMAD::Double &h_max, bool &count_eval)
    {
        char errstr[1024];
        bool stop = false;
        int i, j, n = x.size(), status;
        double *xm, *fvals;
        mxLogical *sur;
        count_eval = true; //mexErrMsgTxt will kill MEX

        //Check for Ctrl-C
        if (utIsInterruptPending()) {
            utSetInterruptPending(false); /* clear Ctrl-C status */
            mexPrintf("\nCtrl-C Detected. Exiting NOMAD...\n\n");
            count_eval = false;
            raise(SIGINT);
            return false;
        }

        //Blackbox / Objective Evaluation
        fun->plhs[0] = NULL;
        xm = mxGetPr(fun->prhs[fun->xrhs]);
        for(i=0;i<n;i++)
            xm[i] = x[i].value();
        //Add Surrogate if present and requested
        if(hasSur) {
            sur=mxGetLogicals(fun->prhs[fun->xrhs+1]);
			(x.get_eval_type()==NOMAD::SGTE)? *sur=true:*sur=false;
        }
        //Call MATLAB Objective
        try {
            status = mexCallMATLAB(1, fun->plhs, fun->nrhs, fun->prhs, fun->f);   
        }
        //Note if these errors occur it is due to errors in MATLAB code, no way to recover?
        catch(exception &e) {
            sprintf(errstr,"Unrecoverable Error from Objective / Blackbox Callback:\n%sExiting NOMAD...\n\n",e.what());
            mexWarnMsgTxt(errstr);
            //Force exit
            raise(SIGINT);
            return false;
        }  
        catch(...) {            
            mexWarnMsgTxt("Unrecoverable Error from Objective / Blackbox Callback, Exiting NOMAD...\n\n");
            //Force exit
            raise(SIGINT);
            return false;
        }
        if(status) {
            mexWarnMsgTxt("Unrecoverable Error from Objective / Blackbox Callback, Exiting NOMAD...\n\n");
            raise(SIGINT);
            return false;
        }
        
        #ifdef OPTI_VERSION  
            //Check we got the correct number back
            if(mxGetM(fun->plhs[0]) != nobj)
                mexErrMsgTxt("Incorrect number of elements returned from the objective function");
            //Set Objective (Or multi-objective)
            fvals = mxGetPr(fun->plhs[0]);
            for(i=0;i<nobj;i++)
                x.set_bb_output(i,fvals[i]);

            //Constraint Evaluation
            if(ncon) {
                //Add Surrogate if present and requested
                if(hasSur) {
                    sur=mxGetLogicals(con->prhs[con->xrhs+1]);
                    (x.get_eval_type()==NOMAD::SGTE)? *sur=true:*sur=false;
                }
                con->plhs[0] = NULL;
                xm = mxGetPr(con->prhs[con->xrhs]);
                for(i=0;i<n;i++)
                    xm[i] = x[i].value();
                //Call MATLAB Constraint
                try {
                    status = mexCallMATLAB(1, con->plhs, con->nrhs, con->prhs, con->f);
                }
                catch(...)
                {
                    mexWarnMsgTxt("Unrecoverable Error from Constraint Callback, Exiting NOMAD...\n\n");
                    //Force exit
                    raise(SIGINT);
                    return false;
                }
                if(status) {
                    mexWarnMsgTxt("Unrecoverable Error from Constraint Callback, Exiting NOMAD...\n\n");
                    raise(SIGINT);
                    return false;
                }
                //Check we got the correct number back
                if(mxGetM(con->plhs[0]) != ncon)
                    mexErrMsgTxt("Incorrect number of elements returned from nonlinear constraint function");
                //Set Constraints
                double *cons = mxGetPr(con->plhs[0]);            
                for(i=0,j=nobj;i<ncon;i++,j++)
                    x.set_bb_output(j,cons[i] - con->nlrhs[i]); //subtract nlrhs        
                // Clean up LHS Ptr
                mxDestroyArray(con->plhs[0]); 
            }     
        #else //GERAD VERSION
            //Check we got the correct number of elements back
            if(mxGetNumberOfElements(fun->plhs[0]) > nobj+ncon)
                mexWarnMsgTxt("Black box returns more elements than required. Please provide a BB_OUTPUT_TYPE consistent with your black box function");
            else if(mxGetNumberOfElements(fun->plhs[0]) < nobj+ncon) {
                mexErrMsgTxt("Insufficient outputs provided by the black box function. Exiting NOMAD...\n\n");
                raise(SIGINT);
                return false;
            }
            //Assign bb output
            fvals = mxGetPr(fun->plhs[0]);
            for(i=0;i<(nobj+ncon);i++)
                x.set_bb_output(i,fvals[i]);
             
        #endif
                    
        //Iteration Callback
        if(iterF->enabled)
        {
            iterF->plhs[0] = NULL;
            memcpy(mxGetData(iterF->prhs[1]), &citer, sizeof(int));
            memcpy(mxGetPr(iterF->prhs[2]), fvals, sizeof(double));
            memcpy(mxGetPr(iterF->prhs[3]), xm, n * sizeof(double));
            try {
                status = mexCallMATLAB(1, iterF->plhs, 4, iterF->prhs, iterF->f);
            }
            catch (...)
            {
                mexWarnMsgTxt("Unrecoverable Error from Iteration Callback, Exiting NOMAD...\n\n");
                //Force exit
                raise(SIGINT);
                return false;
            }  
            if(status) {
                mexWarnMsgTxt("Unrecoverable Error from Iteration Callback, Exiting NOMAD...\n\n");
                raise(SIGINT);
                return false;
            }
            
            //Collect return argument
            stop = *(bool*)mxGetData(iterF->plhs[0]);
            //Clean up Ptr
            mxDestroyArray(iterF->plhs[0]);
        }              
        
        //Add Function Eval Counter
        citer++;
        
        // Clean up LHS Fun Ptr
        mxDestroyArray(fun->plhs[0]);
        
        //Check for iterfun stop
        if(stop)
        {
            mexPrintf("\nIterFun Called Stop. Exiting NOMAD...\n\n");
            count_eval = false;
            raise(SIGINT);
            return false;
        }
        else        
            return true;
    }    
};

//MATLAB MultiObj Evaluator Class
class matlabMEval : public NOMAD::Multi_Obj_Evaluator {
private:
    matlabEval *mEval;
public:
    //Constructor
    matlabMEval(const NOMAD::Parameters &p, usrFcn *_obj, int _nobj, usrFcn *_con, int _ncon, iter_fun_data *_iterF) : NOMAD::Multi_Obj_Evaluator(p)
    {
        mEval = new matlabEval(p,_obj,_nobj,_con,_ncon,_iterF);
    }
    //Destructor
    ~matlabMEval(void)
    {
        delete mEval;
    }
    //Function + Constraint Information
    bool eval_x(NOMAD::Eval_Point &x, const NOMAD::Double &h_max, bool &count_eval)
    {
        return mEval->eval_x(x,h_max,count_eval);
    }        
};

//cout Redirection
struct printfbuf : std::streambuf {
public:
     //Constructor
	 printfbuf() 
	 { 
		 setp(m_buffer, m_buffer + s_size - 2); 
	 }	 
private:
    enum { s_size = 1024 }; //not sure on this size
	char m_buffer[s_size];
	int_type overflow(int_type c) 
	{
		if (!traits_type::eq_int_type(c, traits_type::eof())) {
			*pptr() = traits_type::to_char_type(c);
			pbump(1);
		}
		return sync() != -1 ? traits_type::not_eof(c) : traits_type::eof();
	}

	int sync() {
		*pptr() = 0;
		mexPrintf(pbase());
        mexEvalString("drawnow;");
		setp(m_buffer, m_buffer + s_size - 2);
		return 0;
	} 
 };

// Main Entry Function
// -----------------------------------------------------------------
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    //Input Args
    usrFcn fun, con;
    double *x0, *lb = NULL, *ub = NULL;
    char *xtype = NULL;

    //Outputs Args
    double *x, *fval, *exitflag, *iter, *nfval;
    
    //Internal Vars
    size_t ndec;   
    int i, nobj = 1, ncon = 0;
    char errstr[1024]; //used for returning error info
    iter_fun_data iterF;

    //Check user inputs
    if(!checkInputs(prhs,nrhs,plhs,nlhs))
        return;
    
    //Redirect cout
    printfbuf buf;
    std::streambuf *cout_sbuf = std::cout.rdbuf(); //keep existing buffer
	std::cout.rdbuf(&buf); //redirect buffer
    
    //NOMAD Vars    
    NOMAD::Mads *mads;
    NOMAD::Display out(std::cout);
    NOMAD::Parameters p(out);
    NOMAD::Point px0;
    NOMAD::Double *nx0;
    NOMAD::stop_type stopflag;
    //Evaluator Vars
    matlabEval *mSEval = NULL;
    matlabMEval *mBEval = NULL;
    
    //Set Option Defaults    
    int printLevel = 0;
    char *paramfile = NULL;
    mxArray *bb_out_type = NULL;
    iterF.enabled = false;
 
    //Get Size
    ndec = mxGetNumberOfElements(pX0);
    //Get Blackbox / Objective Function Handle
    if (mxIsChar(pFUN)) {
        if(mxGetString(pFUN, fun.f, FLEN) != 0)
            mexErrMsgTxt("error reading objective name string");
        fun.nrhs = 1;
        fun.xrhs = 0;
    } else {
        fun.prhs[0] = (mxArray*)pFUN;
        strcpy(fun.f, "feval");
        fun.nrhs = 2;
        fun.xrhs = 1;
    }
    fun.prhs[fun.xrhs] = mxCreateDoubleMatrix(ndec, 1, mxREAL); //x
    
    //Get x0
    x0 = mxGetPr(pX0);
    
    //Get xtype
    if(nrhs > eXTYPE && !mxIsEmpty(pXTYPE))
        xtype = mxArrayToString(pXTYPE);
    
    //Get MEX Options if specified
    if(nrhs > eOPTS && !mxIsEmpty(pOPTS)) {
        if(mxGetField(pOPTS,0,"display_degree") && !mxIsEmpty(mxGetField(pOPTS,0,"display_degree")))
            printLevel = (int)*mxGetPr(mxGetField(pOPTS,0,"display_degree"));
        if(mxGetField(pOPTS,0,"param_file") && !mxIsEmpty(mxGetField(pOPTS,0,"param_file")))
            paramfile = mxArrayToString(mxGetField(pOPTS,0,"param_file"));
        if(mxGetField(pOPTS,0,"bb_output_type") && !mxIsEmpty(mxGetField(pOPTS,0,"bb_output_type")))
            bb_out_type = mxGetField(pOPTS,0,"bb_output_type");
        if(mxGetField(pOPTS,0,"iterfun") && !mxIsEmpty(mxGetField(pOPTS,0,"iterfun")))
        {
            iterF.prhs[0] = (mxArray*)mxGetField(pOPTS,0,"iterfun");
            strcpy(iterF.f, "feval");
            iterF.enabled = true;  
            iterF.prhs[1] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
            iterF.prhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
            iterF.prhs[3] = mxCreateDoubleMatrix(ndec,1,mxREAL);
        }
    }     
    
    //Create Outputs (note x and fval created below, due to allowing bi-objective)
    plhs[2] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(1,1, mxREAL);     
    exitflag = mxGetPr(plhs[2]);    
    iter = mxGetPr(plhs[3]);
    nfval = mxGetPr(plhs[4]);
    
    //Setup ndec
    p.set_DIMENSION((int)ndec);  
    //Warn if >1000
    if(ndec > 1000 && printLevel) {
        sprintf(errstr,"Warning: NOMAD is designed problems less than 1000 variables. Your model has %d.\nWhile unlikely, it is possible NOMAD may not perform as intended on this problem.",static_cast<int>(ndec));
        mexWarnMsgTxt(errstr);
    }
    
    //Setup Lower Bounds     
    if(nrhs > eLB && !mxIsEmpty(pLB)) {
        NOMAD::Point ptLB;
        NOMAD::Double *dLB = new NOMAD::Double[ndec];
        lb = mxGetPr(pLB);
        for(i=0;i<ndec;i++) {
            if(!mxIsInf(lb[i])) //if not initialized will not be used
                dLB[i] = lb[i];
        }
        ptLB.set((int)ndec,dLB);
        p.set_LOWER_BOUND(ptLB);
        delete [] dLB;
    }
    //Setup Upper Bounds
    if(nrhs > eUB && !mxIsEmpty(pUB)) {
        NOMAD::Point ptUB;
        NOMAD::Double *dUB = new NOMAD::Double[ndec]; 
        ub = mxGetPr(pUB);
        for(i=0;i<ndec;i++) {
            if(!mxIsInf(ub[i])) //if not initialized will not be used
                dUB[i] = ub[i];
        }
        ptUB.set((int)ndec,dUB);
        p.set_UPPER_BOUND(ptUB);
        delete [] dUB;
    }
    
    //Setup x0
    nx0 = new NOMAD::Double[ndec]; 
    #ifdef OPTI_VERSION
        double xl, xu;
        //If integer variables declared, need to ensure x0[i] is an integer
        if(xtype) {
            for(i=0;i<ndec;i++) {
                switch(tolower(xtype[i]))
                {
                    case 'c':
                        //Ensure within bounds
                        if(lb && x0[i] < lb[i])
                            nx0[i] = lb[i];
                        else if(ub && x0[i] > ub[i])
                            nx0[i] = ub[i];
                        else
                            nx0[i] = x0[i];  //no rounding
                        break;
                    case 'i':
                    case 'b':
                        xl = floor(x0[i]); //First round is a floor
                        //If lower bounds exist
                        if(lb) {
                            //if lower bound broken
                            if(xl < lb[i]) { 
                                xu = ceil(x0[i]);
                                //If upper bounds exist, check bound directions
                                if(ub && xu > ub[i]) { //if broken, no integer x0 exists
                                    sprintf(errstr,"x0[%d] cannot be rounded to an integer value between lb: %g, ub %g",i,lb[i],ub[i]);
                                    mexErrMsgTxt(errstr);
                                }                            
                                if(xu != x0[i]) { //If we changed something, warn user
                                    sprintf(errstr,"x0[%d] was rounded up to %g to suit NOMAD interface",i,xu);
                                    mexWarnMsgTxt(errstr);
                                }
                                //Save ceil value
                                nx0[i] = xu;                            
                            }
                            //Floor value did not break lower bounds, value OK
                            else {
                                if(xl != x0[i]) { //If we changed something, warn user
                                    sprintf(errstr,"x0[%d] was rounded down to %g to suit NOMAD interface",i,xl);
                                    mexWarnMsgTxt(errstr);
                                }
                                //Save floor value
                                nx0[i] = xl;
                            }
                        }
                        //No lower bounds, floor value assumed OK
                        else {
                            if(xl != x0[i]) { //If we changed something, warn user
                                sprintf(errstr,"x0[%d] was rounded down to %g to suit NOMAD interface",i,xl);
                                mexWarnMsgTxt(errstr);
                            }
                            //Save floor value
                            nx0[i] = xl;
                        }
                        break;
                    case 'r':
                        mexErrMsgTxt("Please specify continuous (real) variables using 'c' (as opposed to 'r') when using the OPTI version");
                        break;
                    default:
                        sprintf(errstr,"Unknown xtype[%d] character: %c\n\nValid options are 'C', 'I', or 'B'\n",i,xtype[i]);
                        mexErrMsgTxt(errstr);
                }
            }
        }
        //Else user start position within bounds
        else {
           for(i=0;i<ndec;i++) {
                if(lb && x0[i] < lb[i])
                    nx0[i] = lb[i];
                else if(ub && x0[i] > ub[i])
                    nx0[i] = ub[i];
                else
                    nx0[i] = x0[i]; 
           }
        }
    
    #else //GERAD VERSION - no x0 checking
        for(i=0;i<ndec;i++)
            nx0[i] = x0[i];
    #endif
    px0.set((int)ndec,nx0);
    p.set_X0(px0);
    delete [] nx0;
    
    #ifdef OPTI_VERSION
        //Setup Nonlinear Constraints
        if(nrhs > eNLCON && !mxIsEmpty(pNLCON)) {
            if (mxIsChar(pNLCON)) {
                if(mxGetString(pNLCON, con.f, FLEN) != 0)
                    mexErrMsgTxt("error reading constraint name string");
                con.nrhs = 1;
                con.xrhs = 0;
            } else {
                con.prhs[0] = (mxArray*)pNLCON;
                strcpy(con.f, "feval");
                con.nrhs = 2;
                con.xrhs = 1;
            }
            con.prhs[con.xrhs] = mxCreateDoubleMatrix(ndec, 1, mxREAL); //x
            if(nrhs < eNLRHS+1 || mxIsEmpty(pNLRHS)) {//we will default to <= 0           
                ncon = -1;
                con.nlrhs = NULL;
            }
            else {
                ncon = (int)mxGetNumberOfElements(pNLRHS);       
                con.nlrhs = mxGetPr(pNLRHS);
            }
        }
        //Setup Input Variable Types
        if(xtype)
            p.set_BB_INPUT_TYPE(detInTypes(xtype,ndec));

        //Setup Evaluation Return Types + #'s of Obj + Con
        p.set_BB_OUTPUT_TYPE(detRetTypes(&fun,bb_out_type,&nobj,&con,&ncon,x0,ndec));

        //Set All Normal NOMAD Options
        p.set_DISPLAY_DEGREE(0); //default
        
    #endif //GERAD Version does not have a separate constraint handler and handles input and output types using options
    
    //Set User Options
    if(nrhs > eOPTS && !mxIsEmpty(pOPTS))
        setNOMADopts(p,pOPTS);  
    else
        setNOMADopts(p,NULL);

    //If the user has specified a parameter file to read, see if it exists, and if so, read and parse it.
    if(paramfile) {
        FILE *pFile = fopen(paramfile,"r");
        if(pFile==NULL) {
            sprintf(errstr,"Cannot open parameter file: %s\n\nEnsure it exists!",paramfile);
            mexErrMsgTxt(errstr);
        }
        else{
            fclose(pFile); //close file pointer (we don't need it)
            try
            {
                p.read(paramfile);
            }
            catch(exception &e)
            {
                sprintf(errstr,"NOMAD Parameter File Read Error:\n\n%s",e.what());
                mexErrMsgTxt(errstr);
            }        
        }
    }

    //Check NOMAD parameters
    try {
        p.check();
    }
    catch(exception &e) {        
        sprintf(errstr,"NOMAD Parameter Error:\n\n%s",e.what());
        mexErrMsgTxt(errstr);
    }
    //Check for categorical variables
    if (p.get_signature()->has_categorical())
		mexErrMsgTxt("The MATLAB version of NOMAD does not support categorical variables yet! Check BB_INPUT_TYPE parameter.");

    //If GERAD version, obtain number of objectives and constraints from parameters
    #ifndef OPTI_VERSION
        nobj=p.get_nb_obj();
        ncon=(int)p.get_bb_output_type().size()-nobj;        	
    #endif
        
    //If user has requested surrogates, setup extra input arg to callbacks
    if (p.has_sgte()) {
        fun.prhs[fun.xrhs+1] = mxCreateLogicalMatrix(1,1); //extra logical indicating surrogate or not
        fun.nrhs++;
        #ifdef OPTI_VERSION
            con.prhs[con.xrhs+1] = mxCreateLogicalMatrix(1,1);
            con.nrhs++;
        #endif
    }
        
    //Print Header
    if(printLevel) {
        mexPrintf("\n------------------------------------------------------------------\n");
        mexPrintf(" This is NOMAD v%s\n",NOMAD::VERSION.c_str());
        mexPrintf(" Authors: M. Abramson, C. Audet, G. Couture, J. Dennis,  S. Le Digabel, C. Tribes\n\n");
        mexPrintf(" Problem Properties:\n");
        mexPrintf(" # Decision Variables:               %4d\n",ndec);        
        mexPrintf(" # Number of Objectives:             %4d\n",nobj);
        mexPrintf(" # Number of Nonlinear Constraints:  %4d\n",ncon);

        mexPrintf("------------------------------------------------------------------\n");
        mexEvalString("drawnow;"); //flush draw buffer
    }
    //Let this file sort out errors
    mexSetTrapFlag(1);
    
    //Reset tags and bbe (C.Tribes 3/14)
    NOMAD::Eval_Point::reset_tags_and_bbes();
    
    //Create evaluator and run mads based on number of objectives
    try
    {     
        if(nobj > 1) {
            mBEval = new matlabMEval(p,&fun,nobj,&con,ncon,&iterF); //Bi-Objective Evaluator
            mads = new NOMAD::Mads(p, mBEval); //Run NOMAD  
            stopflag = mads->multi_run();
        }
        else {
            mSEval = new matlabEval(p,&fun,nobj,&con,ncon,&iterF); //Single Objective Evaluator
            mads = new NOMAD::Mads(p, mSEval); //Run NOMAD 
            stopflag = mads->run();
        }
    }
    catch(exception &e)
    {
        //Free Memory (C.Tribes 3/14)
        if(mSEval) delete mSEval; mSEval = NULL;
        if(mBEval) delete mBEval; mBEval = NULL;
        delete mads;
        if(xtype) mxFree(xtype); xtype = NULL;
        //Report Error
        sprintf(errstr,"NOMAD Run Error:\n\n%s",e.what());
        mexErrMsgTxt(errstr);
    }
    
    if(printLevel)
        mexPrintf("------------------------------------------------------------------\n");    
    
    //Obtain Solution based on Single or Multi-Objective (C.Tribes oct 09, 2013 --- changed to deal with pareto output)
    if(nobj>1) {
		NOMAD::Pareto_Front * pareto_front=mads->get_pareto_front();				
		if (pareto_front) {					
			int nb_pareto_pts = pareto_front->size(); 
			plhs[0] = mxCreateDoubleMatrix(ndec,nb_pareto_pts, mxREAL);
			plhs[1] = mxCreateDoubleMatrix(nobj,nb_pareto_pts, mxREAL);  
			x = mxGetPr(plhs[0]); 
			fval = mxGetPr(plhs[1]); 
			const NOMAD::Eval_Point * cur = pareto_front->begin();
			int i=0;
			while ( cur ) {			
                //mexPrintf("i %d OK: %d FEAS: %d\n",i,cur->is_eval_ok(),cur->is_feasible ( p.get_h_min() ));
				if ( cur->is_eval_ok() && cur->is_feasible ( p.get_h_min() ) ) {
					const std::list<int>           & index_obj = p.get_index_obj();
					std::list<int>::const_iterator   it , end  = index_obj.end();
					const NOMAD::Point             & bbo       = cur->get_bb_outputs();
					int                              j         = 0;
					NOMAD::Point multi_obj ( static_cast<int>(index_obj.size()) );
					
					for ( it = index_obj.begin() ; it != end ; ++it,j++ )
						fval[nobj*i+j] = bbo[*it].value();					
					for(j=0;j<ndec;j++)
						x[ndec*i+j] = (*cur)[j].value();
				} 
				cur = pareto_front->next();
				i++;
			}
			*exitflag = getStatus(stopflag);
		}
		else {
			stopflag = (NOMAD::stop_type)10;
			*exitflag = -1; //No solution
		}		
	}
	else {
        //Create single objective outputs
		plhs[0] = mxCreateDoubleMatrix(ndec,1, mxREAL);
		plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL); 
		x = mxGetPr(plhs[0]); 
		fval = mxGetPr(plhs[1]); 		
		const NOMAD::Eval_Point *bestSol = mads->get_best_feasible();
		if(bestSol == NULL) {
			bestSol = mads->get_best_infeasible();  
			//manually set as infeasible (no infeasible stop flag)
			stopflag = (NOMAD::stop_type)10;
		}
		if(bestSol == NULL)        
			*exitflag = -1; //No solution

		//If we have a solution, save it
		if(*exitflag != -1) {
			//Save x
			NOMAD::Point pt(*bestSol);
			for(i=0;i<ndec;i++)
				x[i] = pt[i].value();
			//Save fval
			*fval = bestSol->get_f().value();
			//Save Status
			*exitflag = getStatus(stopflag);			
		}
	}
    //Common outputs
    *iter = mads->get_stats().get_iterations();
	*nfval = mads->get_stats().get_bb_eval();
    
    //Return cout to initial buffer
    std::cout.rdbuf(cout_sbuf);
    //Return error control to default
    mexSetTrapFlag(0);
    
    //Free Memory
    if(mSEval) delete mSEval; mSEval = NULL;
    if(mBEval) delete mBEval; mBEval = NULL;
    delete mads;
    if(xtype) mxFree(xtype); xtype = NULL;
}

//Determine Variable Return Types + # Objectives + #NL Constraints (only OPTI Version)
vector<NOMAD::bb_output_type> detRetTypes(usrFcn *fun, mxArray *out_types, int *nobj, usrFcn *con, int *ncon, double *x0, size_t n)
{
    int i, j, stat;

    //Test Blackbox / Objective Evaluation
    fun->plhs[0] = NULL;
    memcpy(mxGetPr(fun->prhs[fun->xrhs]), x0, n * sizeof(double));
    //Call MATLAB Objective
    stat = mexCallMATLAB(1, fun->plhs, fun->nrhs, fun->prhs, fun->f);
    if(stat)
        mexErrMsgTxt("Error calling Objective Function!");
    //Ensure we have a real column
    if(mxGetN(fun->plhs[0]) > mxGetM(fun->plhs[0]))
        mexErrMsgTxt("The objective function must return a scalar or column vector");
    if(mxIsSparse(fun->plhs[0]) || mxIsComplex(fun->plhs[0]))
         mexErrMsgTxt("The objective function must return a real, dense, vector");
    
    //Ensure 1 or 2 rows (max of bi-objective)
    *nobj = (int)mxGetNumberOfElements(fun->plhs[0]);
    if(*nobj < 1 || *nobj > 2)
        mexErrMsgTxt("The objective function must return a scalar or 2x1 vector (bi-objective)");          
    // Clean up LHS Ptr
    mxDestroyArray(fun->plhs[0]);

    //Test Constraint Evaluation
    if(*ncon) {
        con->plhs[0] = NULL;
        memcpy(mxGetPr(con->prhs[con->xrhs]), x0, n * sizeof(double));
        //Call MATLAB Objective
        stat = mexCallMATLAB(1, con->plhs, con->nrhs, con->prhs, con->f);
        if(stat)
            mexErrMsgTxt("Error calling Constraint Function!");
        //Ensure we have a real column
        if(mxGetN(con->plhs[0]) > mxGetM(con->plhs[0]))
            mexErrMsgTxt("The constraint function must return a scalar or column vector");
        if(mxIsSparse(con->plhs[0]) || mxIsComplex(con->plhs[0]))
            mexErrMsgTxt("The constraint function must return a real, dense, vector");

        //If we have nlrhs, check against returned vector
        if(*ncon > 0) {
            if(mxGetM(con->plhs[0]) != *ncon)
                mexErrMsgTxt("The vector returned from nlcon is not the same length as nlrhs!");
        }
        else {           
            *ncon = (int)mxGetM(con->plhs[0]);        
            con->nlrhs = mxGetPr(mxCreateDoubleMatrix(*ncon,1, mxREAL)); //create dummy rhs
        }
        // Clean up LHS Ptr
        mxDestroyArray(con->plhs[0]);
    }

    //Create Return Vector
    vector<NOMAD::bb_output_type> varType(*nobj+*ncon);
    //Fill Objective Information
    for(i=0;i<*nobj;i++)
        varType[i] = NOMAD::OBJ;
    //If the user has set bb_output_type use it to set constraint types
    if(out_types) {
        if(!mxIsCell(out_types))
            mexErrMsgTxt("Parameter bb_output_type must be a cell array of strings.");
        if(mxGetNumberOfElements(out_types) != *ncon)
            mexErrMsgTxt("You must specify a bb_output_type for each element returned by the constraint function (not objective).");
        //Process each element in the array
        for(i=0,j=*nobj;i<*ncon;i++,j++){
            mxArray *ctype = mxGetCell(out_types,i);
            if(!mxIsChar(ctype))
                mexErrMsgTxt("An element in the bb_output_type paramter is not a string");
            //Get the string
            char *str = mxArrayToString(ctype);
            //Set varType based on string
            varType[j] = getOutEnum(str);
            //Free local memory at each iteration
            mxFree(str);
        } 
    }
    //Or Fill Constraint Information as Default all PB
    else {
        for(i=0,j=*nobj;i<*ncon;i++,j++)
            varType[j] = NOMAD::PB;
    }
  
    return varType;
}

//Convert user string to NOMAD bb_output_type enum
NOMAD::bb_output_type getOutEnum(char *str)
{
    char errstr[1024];
    lower(str); //convert to lowercase
    
    if(!strcmp(str,"obj"))
        #ifdef OPTI_VERSION
            mexErrMsgTxt("The OPTI Version of NOMAD uses bb_output_type to describe CONSTRAINTS ONLY. Please remove any OBJ terms from the cell array.");
        #else
            return NOMAD::OBJ;
        #endif
    else if(!strcmp(str,"pb"))
        return NOMAD::PB;
    else if(!strcmp(str,"eb"))
        return NOMAD::EB;
    else if(!strcmp(str,"peb"))
        return NOMAD::PEB_P;
    else if(!strcmp(str,"f"))
        return NOMAD::FILTER;
    else if(!strcmp(str,"nothing") || !strcmp(str,"-"))
        return NOMAD::UNDEFINED_BBO;
    else {
        sprintf(errstr,"Unknown BB_OUTPUT_TYPE %s",str);
        mexErrMsgTxt(errstr);
    } 
    //Default
    return NOMAD::UNDEFINED_BBO;
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

//Determine input variable types (only OPTI version)
vector<NOMAD::bb_input_type> detInTypes(char *xtype, size_t n)
{
    int i; 
    char msgbuf[1024];
    vector<NOMAD::bb_input_type> varType(n);
    for(i=0;i<n;i++) {
        switch(tolower(xtype[i]))
        {
            case 'c':             
            	varType[i] = NOMAD::CONTINUOUS; break; // C.Tribes C-> categorical (not implemented in Matlab version) R-> continuous/real !
            case 'r':
                mexErrMsgTxt("Please specify continuous (real) variables using 'c' when using the OPTI version");            	
            case 'i':
                varType[i] = NOMAD::INTEGER; break;
            case 'b':
                varType[i] = NOMAD::BINARY; break;
            default:
                sprintf(msgbuf,"Unknown xtype[%d] character: %c\n\nValid options are 'C', 'I', or 'B'\n",i,xtype[i]);
                mexErrMsgTxt(msgbuf);
        }
    }
    return varType;
}

//User Input Checking + Version / Info / Help
int checkInputs(const mxArray *prhs[], int nrhs, mxArray *plhs[], int nlhs)
{    
    size_t ndec;
    char *str = NULL;
    
    //MEX Display Version (OPTI compatibility)
    if (nrhs < 1) {
        if(nlhs < 1)
            printSolverInfo();
        else
            plhs[0] = mxCreateString(NOMAD::VERSION.c_str()); 
        return 0;
    }
    
    //Redirect cout
    printfbuf buf;
    std::streambuf *cout_sbuf = std::cout.rdbuf(); //keep existing buffer
	std::cout.rdbuf(&buf); //redirect buffer
    
    //NOMAD Display
    NOMAD::Display out (std::cout);     
    
    //Check for display on options passed as structure
    if(nrhs == 1 && mxIsStruct(prhs[0])) {     
        int i, no = mxGetNumberOfFields(prhs[0]);
        const char *field;    
        std::list<std::string> ls;    
        //For all fields, display nomad help
        for(i=0;i<no;i++) {
            field = mxGetFieldNameByNumber(prhs[0],i);
            string st(field);
            ls.push_back ( st );
        }
        if(no>0) {
            //NOMAD Display
			NOMAD::Parameters p ( out );
            p.help(ls);            
        }
        std::cout.rdbuf(cout_sbuf); //Return cout to initial buffer
        return 0;
    }
    
    //Check for Version / Information / Help Request
    if(nrhs == 1 && mxIsChar(prhs[0])) {
        str = mxArrayToString(prhs[0]); 
        //Check for Info Request
        if(!strcmp(str,"-I") || !strcmp(str,"-INFO") || !strcmp(str,"-i") || !strcmp(str,"-info")) {
            //Information to Print
            mexPrintf("\nNOMAD Blackbox Optimization Software, v%s\n\n",NOMAD::VERSION.c_str());
            mexPrintf("Copyright (C) 2001-2013\n");
            mexPrintf("Mark A. Abramson     - The Boeing Company\n");
            mexPrintf("Charles Audet        - Ecole Polytechnique de Montreal\n");
            mexPrintf("Gilles Couture       - Ecole Polytechnique de Montreal\n");
            mexPrintf("John E. Dennis, Jr.  - Rice University\n");
            mexPrintf("Sebastien Le Digabel - Ecole Polytechnique de Montreal\n");
			mexPrintf("Christophe Tribes    - Ecole Polytechnique de Montreal\n\n");
            mexPrintf("Funded in part by AFOSR and Exxon Mobil.\n\n");
			mexPrintf("License   : \'%s\'\n",NOMAD::LGPL_FILE.c_str());
			mexPrintf("User guide: \'%s\'\n",NOMAD::USER_GUIDE_FILE.c_str());
			mexPrintf("Examples  : \'%s\'\n",NOMAD::EXAMPLES_DIR.c_str());
			mexPrintf("Tools     : \'%s\'\n\n",NOMAD::TOOLS_DIR.c_str());
			mexPrintf("Please report NOMAD bugs to nomad@gerad.ca and MEX Interface bugs to jocurrie@aut.ac.nz\n\n");
            
            std::cout.rdbuf(cout_sbuf); //Return cout to initial buffer
            return 0;
        }
        //Check for Ver Request
        if(!strcmp(str,"-V") || !strcmp(str,"-v") || !strcmp(str,"-version")) 
		{
            mexPrintf("NOMAD - version %s - www.gerad.ca/nomad\n\n",NOMAD::VERSION.c_str());
            #ifdef OPTI_VERSION
                mexPrintf("MEX Interface (OPTI) v%s\n",NOMADMEX_VERSION);
            #else
                mexPrintf("MEX Interface (GERAD) v%s\n",NOMADMEX_VERSION);
            #endif
            
            std::cout.rdbuf(cout_sbuf); //Return cout to initial buffer            
            return 0;
		}
        //Check for Help Request
        if (strcmp(str,"-H")<0 || strcmp(str,"-HELP")<0 || strcmp(str,"-h")<0 || strcmp(str,"-help")<0 ) 
		{
            NOMAD::Parameters p ( out );
			std::list<std::string> helpS;
			const char * toks=" ";
			char *w = strtok(str,toks) ;
			for ( w = strtok(NULL,toks) ; w != NULL ; w = strtok(NULL,toks) ) 
			{
				helpS.push_back(w);
			}
			p.help (helpS );
            
            std::cout.rdbuf(cout_sbuf); //Return cout to initial buffer
            return 0;
		}
    }
    
    //No NOMAD request, return cout to initial buffer (very Important!)
    std::cout.rdbuf(cout_sbuf); 

    //Otherwise assume we have a normal problem
    if(nrhs < 2)
        mexErrMsgTxt("You must supply at least 2 arguments to nomad!\n\nnomad(fun,x0)\n");
       
    //Check Types
    if(!mxIsFunctionHandle(pFUN) && !mxIsChar(pFUN))
        mexErrMsgTxt("fun must be a function handle or function name!");
    if(!mxIsDouble(pX0) || mxIsComplex(pX0) || mxIsEmpty(pX0))
        mexErrMsgTxt("x0 must be a real double column vector!");
    
    //Get ndec
    ndec = mxGetNumberOfElements(prhs[1]);
    
    //Check Bounds
    if(nrhs > 2) {
        if(!mxIsDouble(pLB) || mxIsComplex(pLB))
            mexErrMsgTxt("lb must be a real double column vector!");
        if(nrhs > 3 && (!mxIsDouble(pUB) || mxIsComplex(pUB)))
            mexErrMsgTxt("ub must be a real double column vector!");
        //Check Sizes
        if(!mxIsEmpty(pLB) && (ndec != mxGetNumberOfElements(pLB)))
            mexErrMsgTxt("lb is not the same length as x0! Ensure they are both Column Vectors");
        if(nrhs > 3 && !mxIsEmpty(pUB) && (ndec != mxGetNumberOfElements(pUB)))
            mexErrMsgTxt("ub is not the same length as x0! Ensure they are both Column Vectors");
    }
    
    #ifdef OPTI_VERSION
        //Version check
        if(nrhs > eNLRHS && mxIsStruct(pNLRHS))    
            mexErrMsgTxt("It appears you may be calling NOMAD using the GERAD NOMAD Syntax - nomad(bb,x0,lb,ub,opts). This is the OPTI NOMAD version and requires nomad(fun,x0,lb,ub,nlcon,nlrhs,xtype,opts).");
    
        //Check Nonlinear Constraint Handle
        if(nrhs > eNLCON && !mxIsEmpty(pNLCON)) {
            if(!mxIsFunctionHandle(pNLCON) && !mxIsChar(pNLCON))
                mexErrMsgTxt("nlcon must be a function handle or function name!");
            if(nrhs > eNLRHS && (!mxIsDouble(pNLRHS) || mxIsComplex(pNLRHS)))
                mexErrMsgTxt("nlrhs must be a real double column vector!");
        }    
        //Check for xtype
        if(nrhs > eXTYPE && !mxIsEmpty(pXTYPE)) {
            if(!mxIsChar(pXTYPE))
                mexErrMsgTxt("The xtype vector must be a char array!");
            if(ndec != mxGetNumberOfElements(pXTYPE))
                mexErrMsgTxt("xtype is not same length as x0! Ensure they are both Column Vectors");
        }
    #else
        //Version check
        if(nrhs > 6)    
            mexErrMsgTxt("It appears you may be calling NOMAD using the OPTI NOMAD Syntax - nomad(fun,x0,lb,ub,nlcon,nlrhs,xtype,opts). This is the GERAD NOMAD version and requires nomad(bb,x0,lb,ub,opts).");  
    #endif              
    
    //Check Options
    if(nrhs > eOPTS && !mxIsEmpty(pOPTS)) {
        //Version check
        if(mxIsDouble(pOPTS))    
            mexErrMsgTxt("It appears you may be calling NOMAD using the OPTI NOMAD Syntax - nomad(fun,x0,lb,ub,nlcon,nlrhs,xtype,opts). This is the GERAD NOMAD version and requires nomad(bb,x0,lb,ub,opts).");
        
        if(!mxIsStruct(pOPTS))
            mexErrMsgTxt("The specified options must be a structure!");
    }
    
    //Return Continue
    return  1;

}

void setNOMADopts(NOMAD::Parameters &p, const mxArray *opts)
{
    char strbuf[1024];
    int i, j, no = 0;
    NOMAD::Parameter_Entries entries;
    const char *field;
    mxArray *value;
    bool doAdd = false;
    size_t m;
    double val;
    bool has_BB_OUTPUT_TYPE = false;
    
    if(opts)
        no = mxGetNumberOfFields(opts);
    
    //For each field, check if it's empty, if not, set it within NOMAD
    for(i=0;i<no;i++) {
        field = mxGetFieldNameByNumber(opts,i);
        value = mxGetFieldByNumber(opts,0,i);
        //Check we don't have an empty or MEX option
        if(!mxIsEmpty(value) && isNMDOption(field)) {
            //Switch on data type
            switch(mxGetClassID(value))
            {
                case mxDOUBLE_CLASS:
                    m = mxGetNumberOfElements(value);
                    //This interface doesn't write vectors as options
                    if(m > 1) {
                        sprintf(strbuf,"This interface does not support vector options for parameter %s. Please provide it as a string",field);
                        mexErrMsgTxt(strbuf);
                    }
                    val = *mxGetPr(value);
                    //Check if we have an integer parameter (I assume this check is ok!)
                    if(floor(val) == val) {
                        sprintf(strbuf,"%s %.0f",field,val); doAdd = true; //write with no decimal points and no exponent
                    }
                    else {//write as double (exponent ok)                    
                        sprintf(strbuf,"%s %g",field,val); doAdd = true;
                    }
                    break;
                    
                case mxCHAR_CLASS:
                    sprintf(strbuf,"%s %s",field,mxArrayToString(value)); 
                    doAdd = true;
                    //GERAD check for BB_OUTPUT_TYPE
                    if(strcmp(field,"BB_OUTPUT_TYPE") == 0 || strcmp(field,"bb_output_type") == 0)
                        has_BB_OUTPUT_TYPE = true;
                    break;                    
                
                case mxCELL_CLASS:
                    #ifdef OPTI_VERSION
                        //Add each string for the specified field
                        m = mxGetNumberOfElements(value);
                        for(j=0;j<m;j++) {
                            mxArray *cellopt = mxGetCell(value,j);
                            if(!mxIsChar(cellopt)) {
                                sprintf(strbuf,"Element %d in the cell array of paramter %s is not a string",j+1,field);
                                mexErrMsgTxt(strbuf);
                            }
                            //Get the string and create the required option string
                            char *str = mxArrayToString(cellopt);
                            sprintf(strbuf,"%s %s",field,str);
                            //Add to our set of parameter entries
                            NOMAD::Parameter_Entry *pe = new NOMAD::Parameter_Entry(strbuf);
                            if(pe->is_ok())
                                entries.insert(pe);
                            else {
                                sprintf(strbuf,"Parameter %s has an error",field);
                                mexErrMsgTxt(strbuf);
                            }
                            //Free local memory at each iteration
                            mxFree(str);
                        }
                    #else
                        sprintf(strbuf,"Parameter %s NOT SET!\nThe GERAD version of the NOMAD interface does not support cell arrays. Please specify them as strings, with spaces between arguments.",field);
                        mexWarnMsgTxt(strbuf);
                    #endif
                    break;
                    
                default:
                    sprintf(strbuf,"Parameter %s NOT SET, CLASSID = %d which is not supported",field,mxGetClassID(value));
                    mexWarnMsgTxt(strbuf);
            }

            //If we have a valid parameter, add it to our set of entries
            if(doAdd) {
                //mexPrintf("NOMAD Option set as: %s\n",strbuf); //enable for debug
                NOMAD::Parameter_Entry *pe = new NOMAD::Parameter_Entry(strbuf);
                if(pe->is_ok())
                    entries.insert(pe);
                else {
                   sprintf(strbuf,"Parameter %s has an error",field);
                   mexErrMsgTxt(strbuf);
                }     
                doAdd = false;
            }
        }
    }
    //GERAD default bb_output_type if not specified
    #ifndef OPTI_VERSION
        if(!has_BB_OUTPUT_TYPE) {
            NOMAD::Parameter_Entry *pe = new NOMAD::Parameter_Entry("BB_OUTPUT_TYPE OBJ");
            entries.insert(pe);
        }
    #endif
    
    try
    {
        p.read(entries);
    }
    catch(exception &e)
    {
        sprintf(strbuf,"NOMAD Parameter Read Error:\n\n%s",e.what());
        mexErrMsgTxt(strbuf);
    }
}

//Separates MEX interface options from NOMAD ones
bool isNMDOption(const char *field)
{
    if(!strcmp(field,"param_file"))
        return false;
    else if(!strcmp(field,"iterfun"))
        return false;
    #ifdef OPTI_VERSION
    else if(!strcmp(field,"bb_output_type"))
        return false;
    #endif
    else
        return true;
}

double getStatus(int stat)
{
    switch((int)stat)
    {        
        case 5:         //mesh minimum
        case 8:         //min mesh size
        case 9:         //min poll size
        case 20:        //ftarget reached
        case 19:        //feas reached
            return 1;
            break;
        case 12:        //max time
        case 13:        //max bb eval
        case 14:        //max sur eval
        case 15:        //max evals
        case 16:        //max sim bb evals
        case 17:        //max iter   
        case 23:        //max multi bb evals
        case 24:        //max mads runs
            return 0;
            break;
        case 10:        //max mesh index
        case 11:        //mesh index limits
        case 18:        //max consec fails
        case 25:        //stagnation
        case 26:        //no pareto
        case 27:        //max cache memory
            return -1;
        case 6:         //x0 fail
        case 7:         //p1_fail            
            return -2;
        case 2:         //Unknown stop reason
            return -3; 
        case 3:         //Ctrl-C
        case 4:         //User-stopped
            return -5;
        default:        //Not assigned flag
            return -3;        
    }
}

//Print Solver Information
void printSolverInfo()
{    
    mexPrintf("\n-----------------------------------------------------------\n");
    #ifdef OPTI_VERSION
        char vbuf[6]; getVSVer(vbuf); 
        mexPrintf(" NOMAD: Nonlinear Optimization using the MADS Algorithm [v%s, Built %s, VS%s]\n",NOMAD::VERSION.c_str(),__DATE__,vbuf);
    #else
        mexPrintf(" NOMAD: Nonlinear Optimization using the MADS Algorithm [v%s, Built %s]\n",NOMAD::VERSION.c_str(),__DATE__);              
    #endif
    mexPrintf("  - Released under the GNU Lesser General Public License: http://www.gnu.org/copyleft/lesser.html\n");
    mexPrintf("  - Source available from: http://www.gerad.ca/NOMAD/Project/Home.html\n");

    mexPrintf("\n MEX Interface J.Currie 2013 (www.i2c2.aut.ac.nz)\n");
    mexPrintf("-----------------------------------------------------------\n");
}