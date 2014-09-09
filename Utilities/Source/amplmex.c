/****************************************************************
Copyright (C) 1997-1998, 2000 Lucent Technologies
All Rights Reserved

Permission to use, copy, modify, and distribute this software and
its documentation for any purpose and without fee is hereby
granted, provided that the above copyright notice appear in all
copies and that both that the copyright notice and this
permission notice and warranty disclaimer appear in supporting
documentation, and that the name of Lucent or any of its entities
not be used in advertising or publicity pertaining to
distribution of the software without specific, written prior
permission.

LUCENT DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.
IN NO EVENT SHALL LUCENT OR ANY OF ITS ENTITIES BE LIABLE FOR ANY
SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER
IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.
****************************************************************/

/* Code is largely modified from the original by J.Currie Oct 2011 
   Copyright (C) 2011-2013 Jonathan Currie (I2C2)
 */

#include "mex.h"
#include "asl_pfgh.h"

//Defines
#define FLEN 512 /* max length of command */
#define pPROB    plhs[0]
#define pSIZE    plhs[1]

#define ASLCMD_ERROR    -1
#define ASLCMD_ISOPEN   0
#define ASLCMD_OPEN     1
#define ASLCMD_CLOSE    2
#define ASLCMD_FUN      5
#define ASLCMD_GRAD     6
#define ASLCMD_CON      7
#define ASLCMD_JAC      8
#define ASLCMD_JACSTR   9
#define ASLCMD_HES      10
#define ASLCMD_HESSTR   11
#define ASLCMD_WRITESOL 14
#define ASLCMD_CONVERT  15

//Structure Arguments
enum{eH,eF,eLB,eUB,eA,eCL,eCU,eQ,eL,eQCIND,eX0,eV0,eSENSE,eOBJBIAS,eCONLIN};

//Macros
#define CHECK(cond, msg) if (!(cond)) { mexErrMsgTxt(msg); }
#define CHECKASL(aslptr) { if(!asl) { mexErrMsgTxt("You have not opened the ASL interface! Use asl('open','file path') first\n"); } }
#define CHECKNRHS(nrhs,no) { if(nrhs < no) { sprintf(msgbuf, "Wrong number of right hand side arguments! Expected %d\n", no); mexErrMsgTxt(msgbuf); } }
//Set error jump (long jumps back to setjmp if ASL detects an error later)
#define SETERRJMP() {nerror = -1; err_jmp = &err_jmp0; what = "(?)"; whatp = &what; if (setjmp(err_jmp0.jb)) { sprintf(msgbuf, "AMPL Solver Library (ASL) had trouble evaluating the %s callback.\n\nThis is normally due to an operation that results in Inf or NaN. Check your initial guess (x0)!\n", *whatp); mexErrMsgTxt(msgbuf); mexExit(); return; } }

//Function Prototypes
void printUtilityInfo();
int getCommand(const mxArray *prhs0);
static void mexExit(void);
static double* sizechk(const mxArray *mp, char *who, int m);
static bool comp_x(double *xbase, double *xnew, int n);

//Global
static char msgbuf[FLEN];            
static double *objx = NULL, *conx = NULL;   //store previous x values used to eval obj and con

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //Possible Inputs
    char fpath[FLEN];
    char msg[FLEN];
    char cmd[FLEN]; //user commmand 
    int sp = 0;
    
    //Outputs
    const char *fnames[15] = {"H","f","lb","ub","A","cl","cu","Q","l","qcind","x0","v0","sense","objbias","conlin"};
    double *sizes;
    
    //Internal Vars
    int ii; size_t i,j,k;       //indexing vars
    char *what, **whatp;        //error message vars
    static FILE *nl;            //file handle
    ASL *asl = cur_ASL;         //Current ASL instance
    int icmd = ASLCMD_ERROR;    //Command Integer
    double *sense;              //Objective sense
    double *objbias;            //Objective bias
    double *con_lin;            //linearity of the constraints (<0 nl, 0 lin, >0 quad)  
    double *isopen;             //Is ASL open
    bool nlcon = false;         //indicates whether any constraint is nonlinear
    double *x;                  //Evaluation point
    double *f, *g, *c = NULL;   //Return pointers
    int nerror;                 //eval errors
    
    //Sparse Indexing
    mwIndex *Ir, *Jc;
    double *Pr;

    //QP Checking Vars
    int nqpz;                   //Number of nz in QP Objective Quadratic Part 
    int nqc_con = 0;            //number of quadratic constraints
    int *QP_ir, *QP_jc;         //Pointers used when calling nqpcheck
    double *QP_pr;
    double *pqi;                //pointer to quadratic index vector
    ograd *og;                  //objective gradient structure
            
    //Jacobian Vars
    static double *J = NULL;        //Memory to store intermediate Jacobian Values when using Dense Mode
    static double *J1 = NULL;       //Memory to store Jacobian Values 
    cgrad *cg, **cgp, **cgpe;       //constraint gradient structures
    int *cs;                        //Column starts
    
    //Hessian Vars
    static double *Hsp = NULL;      //Memory to store Hessian Values
    static int nhnz;                //Number of Hessian nz
    double *s, *v;                  //Sigma, Lambda
	int *hcs, *hr;                  //Hessian column starts, row indexs
	double *H, *He,  *W;    	      
    
    //Error catching
    Jmp_buf err_jmp0;
    
    //If no inputs, just return info
    if(nrhs < 1) {
        printUtilityInfo();
        return;
    }
        
    //Get User Command
    icmd = getCommand(prhs[0]);
    
    //Switch Yard for Command
    switch(icmd)
    {
        case ASLCMD_ISOPEN:
            isopen = mxGetPr(plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL));
            if(asl)
                *isopen = 1;
            else
                *isopen = 0;
            break;
        
        case ASLCMD_OPEN:
            //Check for Errors
            if(nrhs < 2)
                mexErrMsgTxt("Expected two arguments to open a file! [x0,v0,lb,ub,cl,cu,sense,sizes] = asl('open','file path')\n");
            if(!mxIsChar(prhs[1]))
                mexErrMsgTxt("File path must be a char array!");
            //Get String
            CHECK(mxGetString(prhs[1], fpath, FLEN) == 0,"error reading file path!");
            //Clear any existing objects
            if (cur_ASL)
                ASL_free(&cur_ASL);
            //Set MEX exit function
            mexAtExit(mexExit);
            
            //Open file for LP/QP/QCQP checking
            asl = ASL_alloc(ASL_read_fg);               //allocate for qp read
            return_nofile = 1;                          //return 0 if stub doesn't exist
            nl = jac0dim(fpath,(ftnlen)strlen(fpath));  //read in passed file
            //Check we got the file
            if(!nl) {
                sprintf(msgbuf, "Can't open (or error opening) %s\n", fpath);
                mexErrMsgTxt(msgbuf);
			}
            //Allocate Vector Memory
            pPROB = mxCreateStructMatrix(1,1,15,fnames);
            mxSetField(pPROB,0,fnames[eX0],mxCreateDoubleMatrix(n_var,1, mxREAL));      
            mxSetField(pPROB,0,fnames[eV0],mxCreateDoubleMatrix(n_con, 1, mxREAL));
            mxSetField(pPROB,0,fnames[eLB],mxCreateDoubleMatrix(n_var, 1, mxREAL));
            mxSetField(pPROB,0,fnames[eUB],mxCreateDoubleMatrix(n_var, 1, mxREAL));            
            mxSetField(pPROB,0,fnames[eCL],mxCreateDoubleMatrix(n_con, 1, mxREAL));
            mxSetField(pPROB,0,fnames[eCU],mxCreateDoubleMatrix(n_con, 1, mxREAL));
            mxSetField(pPROB,0,fnames[eSENSE],mxCreateDoubleMatrix(1, 1, mxREAL));
            mxSetField(pPROB,0,fnames[eOBJBIAS],mxCreateDoubleMatrix(1, 1, mxREAL));
            mxSetField(pPROB,0,fnames[eCONLIN],mxCreateDoubleMatrix(n_con, 1, mxREAL));
            //Get Fields (ASL will fill)       
            X0 = mxGetPr(mxGetField(pPROB,0,fnames[eX0]));   
            pi0 = mxGetPr(mxGetField(pPROB,0,fnames[eV0]));  
            LUv = mxGetPr(mxGetField(pPROB,0,fnames[eLB]));  
            Uvx = mxGetPr(mxGetField(pPROB,0,fnames[eUB]));              
            LUrhs = mxGetPr(mxGetField(pPROB,0,fnames[eCL]));  
            Urhsx = mxGetPr(mxGetField(pPROB,0,fnames[eCU]));  
            sense = mxGetPr(mxGetField(pPROB,0,fnames[eSENSE])); 
            objbias = mxGetPr(mxGetField(pPROB,0,fnames[eOBJBIAS]));
            con_lin = mxGetPr(mxGetField(pPROB,0,fnames[eCONLIN]));  
            //Other Output Args
            sizes = mxGetPr(pSIZE = mxCreateDoubleMatrix(16, 1, mxREAL));
                     
            //Check for complementarity problems
            if(n_cc)
                mexWarnMsgTxt("Ignoring Complementarity Constraints!");
            //Assign asl problem sizes
            sizes[0] = (double)n_var; sizes[1] = (double)n_con; sizes[2] = (double)nzc;
            sizes[3] = (double)lnc; sizes[4] = (double)nbv; sizes[5] = (double)niv;
            sizes[6] = (double)nlc; sizes[7] = (double)nlnc; sizes[8] = (double)nlo;
            sizes[9] = (double)nlvb; sizes[10] = (double)nlvc; sizes[11] = (double)nlvo;
            sizes[12] = (double)nlvbi; sizes[13] = (double)nlvci; sizes[14] = (double)nlvoi;
            sizes[15] = (double)nwv; 
            //Read In For QP Checking
            qp_read(nl,0); 
            //Assign sense
            if(objtype[0] == 1)
                *sense = -1; //max
            else
                *sense = 1; //min            
            
            //Determine Objective Linearity
            nqpz = nqpcheck(0, &QP_ir, &QP_jc, &QP_pr); //check objective for qp
            //Determine Constraints Linearity
            for(ii = 0; ii < n_con; ii++) {
                con_lin[ii] = (double)nqpcheck(-(ii+1), &QP_ir, &QP_jc, &QP_pr);
                if(con_lin[ii] < 0)
                    nlcon = true;
                else if(con_lin[ii] > 0)
                {
                    //nqpz indicates quadratic constraint, ensure is inequality
                    if(LUrhs[ii] != Urhsx[ii])
                        nqc_con++;
                    else
                        nlcon = true; //quadratic equalities not currently handled by any explicit QCQP solver (I know of), make nl
                }                    
            }
            //Check to force for NL
            if(nrhs > 2 && *mxGetPr(prhs[2])==1)
                nlcon = true;
            
            //If objective or any constraint is nonlinear, then we have to process as an NLP
            if(nqpz < 0 || nlcon) {
                //Free the QP read memory
                ASL_free(&asl);
                //Re-open for full NLP read
                asl = ASL_alloc(ASL_read_pfgh);                 //allocate memory for pfgh read
                nl = jac0dim(fpath,(ftnlen)strlen(fpath));      //read passed file (full nl read)
                //Allocate Jacobian Memory [note use M1alloc to let ASL clean it up if multiple instances opened]
                J = (double*)M1alloc(nzc*sizeof(double));       //Memory to store Jacobian nzs  
                //Assign memory for saving obj + con x
                objx = (double*)M1alloc(n_var*sizeof(double));
                conx = (double*)M1alloc(n_var*sizeof(double));
                //Read File (f + g + H)
                pfgh_read(nl, ASL_findgroups); 
                //Assign Hessian Memory
                nhnz = sphsetup(1, 1, n_con > 0, 0);            //one obj, use sigma, optionally use lambda, full hessian
                Hsp = (double*)M1alloc(nhnz*sizeof(double));    //memory to store hessian nzs
            }
            //Otherwise we can process as a LP, QP or QCQP
            else {                
                //Re open for new QP read [appears we can only use nqpcheck once]
                nl = jac0dim(fpath,(ftnlen)strlen(fpath));  //read in passed file
                //Assign objective bias
                *objbias = objconst(0);
                qp_read(nl,0); 
                //Check for quadratic objective
                if(nqpz > 0) {
                    //Capture Pointers
                    nqpz = nqpcheck(0, &QP_ir, &QP_jc, &QP_pr); //check objective for qp
                    //Create QP H
                    mxSetField(pPROB,0,fnames[eH],mxCreateSparse(n_var,n_var,nqpz,mxREAL));                   
                    //Copy in Objective Quadratic Elements (copy-cast where appropriate)
                    memcpy(mxGetPr(mxGetField(pPROB,0,fnames[eH])),QP_pr,nqpz*sizeof(double));
                    Jc = mxGetJc(mxGetField(pPROB,0,fnames[eH]));
                    Ir = mxGetIr(mxGetField(pPROB,0,fnames[eH]));
                    for(i = 0; i <= n_var; i++)
                        Jc[i] = (mwIndex)QP_jc[i];
                    for(i = 0; i < nqpz; i++)
                        Ir[i] = (mwIndex)QP_ir[i];                       
                }
                else //create an empty sparse matrix
                    mxSetField(pPROB,0,fnames[eH],mxCreateSparse(n_var,n_var,0,mxREAL));
                
                //Create QP f
                mxSetField(pPROB,0,fnames[eF],mxCreateDoubleMatrix(n_var,1,mxREAL));
                Pr = mxGetPr(mxGetField(pPROB,0,fnames[eF]));
                //Copy in Objective Linear Elements
                for( og = Ograd[0]; og; og = og->next )
                    Pr[og->varno] = og->coef;
                
                //Create A (linear constraints)
                mxSetField(pPROB,0,fnames[eA],mxCreateSparse(n_con, n_var, nzc, mxREAL));
                if(n_con) {
                    Pr = mxGetPr(mxGetField(pPROB,0,fnames[eA]));
                    Ir = mxGetIr(mxGetField(pPROB,0,fnames[eA]));;                    
                    //Fill in A (will double on quadratic linear sections, but easier to remove once in MATLAB)
                    for(Jc = mxGetJc(mxGetField(pPROB,0,fnames[eA])), cs = A_colstarts, i = 0; i <= n_var; ++i)
                        Jc[i] = (mwIndex)cs[i];
                    cgp = Cgrad;
                    for(i = 0; i < n_con; i++)
                        for(cg = *cgp++; cg; cg = cg->next) {
                            Ir[cg->goff] = (mwIndex)i; 
                            Pr[cg->goff] = cg->coef;
                        }
                }
                
                //Add quadratic constraints if present
                if(nqc_con) {
                    //Allocate a Cell Array to store the quadratic constraint Qs, and vector to store indices
                    mxSetField(pPROB,0,fnames[eQ],mxCreateCellMatrix(nqc_con,1)); //Q
                    mxSetField(pPROB,0,fnames[eL],mxCreateDoubleMatrix(n_var, nqc_con,mxREAL)); //l
                    mxSetField(pPROB,0,fnames[eQCIND],mxCreateDoubleMatrix(nqc_con,1,mxREAL)); //ind                   
                    pqi = mxGetPr(mxGetField(pPROB,0,fnames[eQCIND]));
                    //Fill In Constraints
                    for(ii=0,j=0;ii<n_con;ii++) {
                        //Quadratic Constraints
                        if(con_lin[ii] > 0) {
                            //Create index
                            pqi[j] = ii+1; //increment for matlab index
                            //Capture Pointers
                            nqpz = nqpcheck(-(ii+1), &QP_ir, &QP_jc, &QP_pr); //check constraint for qp;
                            if(nqpz <= 0)
                                mexErrMsgTxt("Error reading quadratic constraints. Assumed constraint was quadratic based on prescan, now appears not?");
                            //Create QC Q
                            mxSetCell(mxGetField(pPROB,0,fnames[eQ]),j,mxCreateSparse(n_var,n_var,nqpz,mxREAL));                   
                            //Copy in Constraint Quadratic Elements (copy-cast where appropriate)
                            Pr = mxGetPr(mxGetCell(mxGetField(pPROB,0,fnames[eQ]),j));
                            Jc = mxGetJc(mxGetCell(mxGetField(pPROB,0,fnames[eQ]),j));
                            Ir = mxGetIr(mxGetCell(mxGetField(pPROB,0,fnames[eQ]),j));
                            for(k = 0; k <= n_var; k++)
                                Jc[k] = (mwIndex)QP_jc[k];
                            for(k = 0; k < nqpz; k++) {
                                Ir[k] = (mwIndex)QP_ir[k];
                                Pr[k] = 0.5*QP_pr[k];  //to QP form
                            }
                            //Create QC l (not sure why we can't extract this from Jacobian, values are wrong)
                            Pr = mxGetPr(mxGetField(pPROB,0,fnames[eL]));
                            for( cg = Cgrad[ii]; cg; cg = cg->next )
                                Pr[j*n_var + cg->varno] = cg->coef;
                            //Increment for next cell / col
                            j++;
                        }
                    } 
                }
                //Put back into function eval mode (just in case)
                qp_opify();
                
            }
            break;
            
        case ASLCMD_CLOSE:
            //Check for Errors
            CHECKASL(asl);
            //Call Exit Function
            mexExit();          
            break;                    
            
        case ASLCMD_FUN:
            //Check for Errors
            CHECKASL(asl);
            CHECKNRHS(nrhs,2);             
            //Get x and check dimensions
            x = sizechk(prhs[1],"x",n_var); 
            //Save x
            if(objx) memcpy(objx,x,n_var*sizeof(double));                   
            //Create objective val memory and get it from ASL       
            SETERRJMP(); what = "objective";            
			f = mxGetPr(plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL));            
			*f = objval(0, x, &nerror);        
            break;
            
        case ASLCMD_GRAD:
            //Check for Errors
            CHECKASL(asl);
            CHECKNRHS(nrhs,2);            
            //Get x and check dimensions
            x = sizechk(prhs[1],"x",n_var);
            //Save x
            if(objx) memcpy(objx,x,n_var*sizeof(double));            
            //Create objective grad memory and get it from ASL     
            SETERRJMP(); what = "gradient";            
			g = mxGetPr(plhs[0] = mxCreateDoubleMatrix(1, n_var, mxREAL));            
			objgrd(0, x, g, &nerror);            
            break;
            
        case ASLCMD_CON:
            //Check for Errors
            CHECKASL(asl);
            CHECKNRHS(nrhs,2);            
            //Get x and check dimensions
            x = sizechk(prhs[1],"x",n_var);
            //Save x
            if(conx) memcpy(conx,x,n_var*sizeof(double));                        
            //Create constraint memory and get it from ASL  
            SETERRJMP(); what = "constraints";
			c = mxGetPr(plhs[0] = mxCreateDoubleMatrix(n_con, 1, mxREAL));   
            if(n_con)
                conval(x, c, &nerror);            
            break;
            
        case ASLCMD_JAC:
            //Check for Errors
            CHECKASL(asl);
            CHECKNRHS(nrhs,2);   
            //Get x and check dimensions
            x = sizechk(prhs[1],"x",n_var);
            //Save x
            if(conx) memcpy(conx,x,n_var*sizeof(double));            
            //Create constraint jac memory and get it from ASL
            SETERRJMP(); what = "Jacobian";            
            //Check for sparsity
            if(nrhs > 2 && *mxGetPr(prhs[2])) {
                sp = 1;
                J1 = mxGetPr(plhs[0] = mxCreateSparse(n_con, n_var, nzc, mxREAL));
            }
            else {
                sp = 0;
                J1 = mxGetPr(plhs[0] = mxCreateDoubleMatrix(n_con, n_var, mxREAL));
            }        
            //Evaluate if we have constraints
            if (n_con) {                
                //Sparse
                if(sp) {
                    jacval(x, J1, &nerror);
                    Ir = mxGetIr(plhs[0]);
                    for(Jc = mxGetJc(plhs[0]), cs = A_colstarts, i = 0; i <= n_var; ++i)
                        Jc[i] = (mwIndex)cs[i];
                    cgp = Cgrad;
                    for(i = 0; i < n_con; i++)
                        for(cg = *cgp++; cg; cg = cg->next)
                            Ir[cg->goff] = (mwIndex)i;  
                }
                //Dense
                else {      
                    jacval(x, J, &nerror);
                    cgp = Cgrad;
                    for(cgpe = cgp + n_con; cgp < cgpe; J1++)
                        for(cg = *cgp++; cg; cg = cg->next)
                            J1[n_con*cg->varno] = J[cg->goff];
                }
            }                        
            break;
            
        case ASLCMD_JACSTR:
            //Check for Errors
            CHECKASL(asl);
            CHECKNRHS(nrhs,1);            
            //Create constraint jacstr memory and get it from ASL
            SETERRJMP(); what = "Jacobian Structure)";                       
            J1 = mxGetPr(plhs[0] = mxCreateSparse(n_con, n_var, nzc, mxREAL));
            //Fill In Structure
            for(i=0;i<nzc;i++)
                J1[i] = 1.0;
            for(Jc = mxGetJc(plhs[0]), cs = A_colstarts, i = 0; i <= n_var; ++i)
                Jc[i] = (mwIndex)cs[i];
            cgp = Cgrad;
            Ir = mxGetIr(plhs[0]);
            for(i = 0; i < n_con; i++)
                for(cg = *cgp++; cg; cg = cg->next)
                    Ir[cg->goff] = (mwIndex)i;                
            break;
            
        case ASLCMD_HES:
            //Check for Errors
            CHECKASL(asl);
            CHECKNRHS(nrhs,4); //assume hess(x,sigma,lambda) and optionally sparse            
            //Check dimensions & get args
            x = sizechk(prhs[1],"x",n_var);
            s = sizechk(prhs[2],"sigma",1);
            v = sizechk(prhs[3],"lambda",n_con);
            
            //Check for sparsity
            if(nrhs > 4 && *mxGetPr(prhs[4])) {
                sp = 1;
                W = mxGetPr(plhs[0] = mxCreateSparse(n_var, n_var, nhnz, mxREAL));
            }
            else {
                sp = 0;    
                W = mxGetPr(plhs[0] = mxCreateDoubleMatrix(n_var, n_var, mxREAL));
            }
            //Check if we need to recalculate objective / constraints
            if(!comp_x(objx,x,n_var)) {
                //Setup Error Catching
                SETERRJMP(); what = "Objective for Hessian";                
                //Re-evaluate Objective
                objval(0, x, &nerror);                
            }            
            if(!comp_x(conx,x,n_var)){
                if(!c)
                    c = mxGetPr(mxCreateDoubleMatrix(n_con, 1, mxREAL));                
                //Setup Error Catching
                SETERRJMP(); what = "Constraints for Hessian";                
                //Re-evaluate Constraints
                conval(x, c, &nerror);
            }            
            //Setup Error Catching
            SETERRJMP(); what = "Hessian";
            
            //Sparse
            if(sp) {
                //This function returns the full (symmetric) Hessian as setup above
                sphes(H = Hsp, 1, s, v);                
                Ir = mxGetIr(plhs[0]);
                Jc = mxGetJc(plhs[0]);
                hcs = sputinfo->hcolstarts;
                hr = sputinfo->hrownos;
                for(i = 0; i <= n_var; i++)
                    Jc[i] = (mwIndex)hcs[i];
                He = H + hcs[n_var];
                while(H < He) {
                    *W++ = *H++;
                    *Ir++ = (mwIndex)*hr++;
                }	
            }
            //Dense
            else
                fullhes(W, n_var, 1, s, v);            
            break;
            
        case ASLCMD_HESSTR:
            //mexPrintf("CMD: Get Hessian Structure\n");
            //Check for Errors
            CHECKASL(asl);
            CHECKNRHS(nrhs,1);            
            //Create hessianstr memory and get it from ASL
            SETERRJMP(); what = "Hessian Structure";
            W = mxGetPr(plhs[0] = mxCreateSparse(n_var, n_var, nhnz, mxREAL));
            Ir = mxGetIr(plhs[0]);
            Jc = mxGetJc(plhs[0]);
            //Get Sparse Info
            hcs = sputinfo->hcolstarts;
            hr = sputinfo->hrownos;
            //Assign col starts
            for(i = 0; i <= n_var; i++)
                Jc[i] = (mwIndex)hcs[i];
            //Assign rows + 1.0 for nz positions
            H = Hsp;                //Start of nz Hsp elements
            He = H + hcs[n_var];    //End of nz Hsp elements
            while(H < He) {
                *W++ = 1.0;                
                *Ir++ = (mwIndex)*hr++;
                *H++; //increment nz element position
            }	                        
            break;           
            
        case ASLCMD_WRITESOL:
            //Check for Errors
            CHECKASL(asl);
            CHECKNRHS(nrhs,2); //asl('writesol',msg,x)            
            //Get Input Args
            CHECK(mxGetString(prhs[1], msg, FLEN) == 0,"error reading message!");
            x = sizechk(prhs[2],"x",n_var);            
            //Write to solution stub file
            write_sol(msg,x,NULL,NULL);
            break;
            
        default:
            mexExit(); //clean up
            mxGetString(prhs[0], cmd, FLEN);
            sprintf(msgbuf, "ASL Command Error! Unknown Command: '%s'\n", cmd);
            mexErrMsgTxt(msgbuf);
            break;
    }
}

//Read user command string and return integer define
int getCommand(const mxArray *prhs0)
{
    char cmd[FLEN]; //user commmand 
    
    if(!mxIsChar(prhs0))
        mexErrMsgTxt("Command must be a lowercase char array!");
    
    CHECK(mxGetString(prhs0, cmd, FLEN) == 0,"error reading command!");
    
    if(!strcmp(cmd,"open"))
        return ASLCMD_OPEN;
    else if(!strcmp(cmd,"isopen"))
        return ASLCMD_ISOPEN;
    else if(!strcmp(cmd,"close"))
        return ASLCMD_CLOSE;
    else if(!strcmp(cmd,"fun") || !strcmp(cmd,"obj"))
        return ASLCMD_FUN;
    else if(!strcmp(cmd,"grad"))
        return ASLCMD_GRAD;
    else if(!strcmp(cmd,"con"))
        return ASLCMD_CON;
    else if(!strcmp(cmd,"jac"))
        return ASLCMD_JAC;
    else if(!strcmp(cmd,"jacstr"))
        return ASLCMD_JACSTR;
    else if(!strcmp(cmd,"hes") || !strcmp(cmd,"hess"))
        return ASLCMD_HES;
    else if(!strcmp(cmd,"hesstr") || !strcmp(cmd,"hessstr"))
        return ASLCMD_HESSTR;
    else if(!strcmp(cmd,"writesol") || !strcmp(cmd,"write_sol"))
        return ASLCMD_WRITESOL;
    else if(!strcmp(cmd,"convert"))
        return ASLCMD_CONVERT;
    else
        return ASLCMD_ERROR;
}

//Called when MEX file is cleared or user exits Matlab (or close is called)
static void mexExit(void)
{
    if (cur_ASL)
        ASL_free(&cur_ASL);    
}

static double* sizechk(const mxArray *mp, char *who, int m)
{
    size_t nel;
    nel = mxGetNumberOfElements(mp); 
    if (nel != m) {
        sprintf(msgbuf,"Expected %s to be %d x 1 rather than %d x %d\n", who, m, mxGetM(mp), mxGetN(mp));
        mexErrMsgTxt(msgbuf);
    }
    return mxGetPr(mp);
}

//Compare x to determine if new obj / con eval required for Hessian
static bool comp_x(double *xbase, double *xnew, int n)
{
    int i;
    //Loop and compare each (should be identical, no need for tolerance)
    for(i = 0; i < n; i++)
    {
        if(xbase[i] != xnew[i])
            return false;
    }
    //Default same
    return true;    
}

void printUtilityInfo()
{
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" ASL: AMPL Solver Library [v%d, Built %s]\n",ASLdate_ASL,__DATE__);
    mexPrintf("  - Source available from: http://www.netlib.org/ampl/solvers/\n");

    mexPrintf("\n MEX Interface J.Currie 2013 (www.i2c2.aut.ac.nz)\n");
    mexPrintf("-----------------------------------------------------------\n");
}