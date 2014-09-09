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
   Copyright (C) 2011 Jonathan Currie (I2C2)
 */

/* Note I am not 100% sure how to deal with mxCreate() functions and their
 * memory requirements when the mex routine does not close. Anyone with
 * experience here and could take a look would be appreciated! */

#include "mex.h"
#include "asl_pfgh.h"

//Defines
#define FLEN 512 /* max length of command */

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
#define ASLCMD_QPHES    12

#define ASLCMD_WRITESOL 14
#define ASLCMD_CONVERT  15

//Macros
#define CHECK(cond, msg) if (!(cond)) { mexErrMsgTxt(msg); }
#define CHECKASL(aslptr) { if(!asl) { mexErrMsgTxt("You have not opened the ASL interface! Use asl('open','file path') first\n"); } }
#define CHECKNRHS(nrhs,no) { if(nrhs < no) { sprintf(msgbuf, "Wrong number of right hand side arguments! Expected %d\n", no); mexErrMsgTxt(msgbuf); } }
//Set error jump (long jumps back to setjmp if ASL detects an errror later)
#define SETERRJMP() {nerror = -1; err_jmp = &err_jmp0; what = "(?)"; whatp = &what; if (setjmp(err_jmp0.jb)) { sprintf(msgbuf, "AMPL Solver Library (ASL) had trouble evaluating the %s callback.\n\nThis is normally due to an operation that results in Inf or NaN. Check your initial guess (x0)!\n", *whatp); mexErrMsgTxt(msgbuf); mexExit(); return; } }

//Function Prototypes
void printUtilityInfo();
int getCommand(const mxArray *prhs0);
static void mexExit(void);
static double* sizechk(const mxArray *mp, char *who, fint m);
static bool comp_x(double *xbase, double *xnew, fint n);

//Global
static char msgbuf[FLEN];
            
static mwIndex *jc = NULL, *ir = NULL;
static double *pr = NULL;

static double *objx = NULL, *conx = NULL;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //Possible Inputs
    char fpath[FLEN];
    char msg[FLEN];
    char cmd[FLEN]; //user commmand 
    
    //Internal Vars
    char *what, **whatp;        //error message vars
    static FILE *nl;            //file handle
    ASL *asl = cur_ASL;         //Current ASL instance
    int icmd = ASLCMD_ERROR;    //Command Integer
    static fint n, nc, nhnz, nz, nelq; //Problem Sizes
    double *isopen;
    
    //Jac + Hess Stuff 
    static double *J, *Hsp;      
    static size_t Jsize;
    size_t i;
    mwIndex *Ir, *Jc;    
	fint *hcs, *hr, nerror;
	double *H, *He, *J1, *c = NULL, *f, *g, *x, *W, *v, *s;
    double *sense, *sizes;
	cgrad *cg, **cgp, **cgpe;
    int sp = 0, *cs;
    
    //QP stuff
    static fint nqpz = 0;
    static double *delsqp; 
    static fint *rowqp, *colqp;
    static mxArray *qpH;
        	
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
            //mexPrintf("CMD: Open File\n");
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
                            
            //Init QP stuff
            jc = NULL; ir = NULL; pr = NULL; 
            
            //Do (QC)QP Check First
            asl = ASL_alloc(ASL_read_fg);               //allocate for qp read
            return_nofile = 1;                          //return 0 if stub doesn't exist
            nl = jac0dim(fpath,(ftnlen)strlen(fpath));  //read in passed file
            //Check we got the file
            if(!nl) {
                sprintf(msgbuf, "Can't open (or error opening) %s\n", fpath);
                mexErrMsgTxt(msgbuf);
			}
            X0 = (double*)M1alloc(n_var*sizeof(double));
            qp_read(nl,0); //qp read
            nqpz = nqpcheck(0, &rowqp, &colqp, &delsqp); //check objective for qp
            //Save QP Hessian if we found one
            if(nqpz > 0) {
                //Copy QP Stuff (create a static copy as colqp etc will be cleaned up)
                jc = malloc((n_var+1)*sizeof(mwIndex)); 
                ir = malloc(nqpz*sizeof(mwIndex));
                pr = malloc(nqpz*sizeof(double));
                for(i = 0; i <= n_var; i++)
                    jc[i] = (mwIndex)colqp[i];
                for(i = 0; i < nqpz; i++)
                    ir[i] = (mwIndex)rowqp[i];                
                memcpy(pr,delsqp,nqpz*sizeof(double));
            }            
            mexPrintf("ncon: %d\n",n_con);
            ASL_free(&asl);
            
            //Attempt to open with ASL
            asl = ASL_alloc(ASL_read_pfgh);         //allocate memory for pfgh read
            nl = jac0dim(fpath,(ftnlen)strlen(fpath));      //read passed file (full nl read)       
            //Read in sizes
            n = n_var;
            nc = n_con;
            nz = nzc;  
            //Check for complementarity problems
            if(n_cc)
                mexWarnMsgTxt("Ignoring Complementarity Constraints!");
            //Assign Jac Memory
            J = (double*)M1alloc(nz*sizeof(double));
            Jsize = nc*n*sizeof(double);
            //Get Memory for Other Args & assign to ASL
            X0 = mxGetPr(plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL));
            pi0 = mxGetPr(plhs[1] = mxCreateDoubleMatrix(nc, 1, mxREAL));
            LUv = mxGetPr(plhs[2] = mxCreateDoubleMatrix(n, 1, mxREAL));
            Uvx = mxGetPr(plhs[3] = mxCreateDoubleMatrix(n, 1, mxREAL));            
            LUrhs = mxGetPr(plhs[4] = mxCreateDoubleMatrix(nc, 1, mxREAL));
            Urhsx = mxGetPr(plhs[5] = mxCreateDoubleMatrix(nc, 1, mxREAL));
            sense = mxGetPr(plhs[6] = mxCreateDoubleMatrix(1, 1, mxREAL));
            sizes = mxGetPr(plhs[7] = mxCreateDoubleMatrix(17, 1, mxREAL));
            //Assign memory for obj + con x
            objx = (double*)malloc(n*sizeof(double));
            conx = (double*)malloc(n*sizeof(double));
            //Read File
            pfgh_read(nl, ASL_findgroups); 
            //Assign sense
            if(objtype[0] == 1)
                *sense = -1; //max
            else
                *sense = 1; //min
            //Assign sizes
            sizes[0] = (double)n_var; sizes[1] = (double)n_con; sizes[2] = (double)nzc;
            sizes[3] = (double)lnc; sizes[4] = (double)nbv; sizes[5] = (double)niv;
            sizes[6] = (double)nlc; sizes[7] = (double)nlnc; sizes[8] = (double)nlo;
            sizes[9] = (double)nlvb; sizes[10] = (double)nlvc; sizes[11] = (double)nlvo;
            sizes[12] = (double)nlvbi; sizes[13] = (double)nlvci; sizes[14] = (double)nlvoi;
            sizes[15] = (double)nwv; sizes[16] = (double)nqpz;

            //Assign Hessian Memory
            nhnz = sphsetup(1, 1, nc > 0, 0);
            Hsp = (double*)M1alloc(nhnz*sizeof(double)); 
            break;
            
        case ASLCMD_CLOSE:
            //mexPrintf("CMD: Close File\n");
            //Check for Errors
            CHECKASL(asl);
            //Call Exit Function
            mexExit();          
            break;                    
            
        case ASLCMD_FUN:
            //mexPrintf("CMD: Get Function Val\n");
            //Check for Errors
            CHECKASL(asl);
            CHECKNRHS(nrhs,2); 
            
            //Get x and check dimensions
            x = sizechk(prhs[1],"x",n); 
            //Save x
            memcpy(objx,x,n*sizeof(double));                   
            //Create objective val memory and get it from ASL       
            SETERRJMP();
			f = mxGetPr(plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL));
            what = "objective";
			*f = objval(0, x, &nerror);        
            break;
            
        case ASLCMD_GRAD:
            //mexPrintf("CMD: Get Gradient Val\n");
            //Check for Errors
            CHECKASL(asl);
            CHECKNRHS(nrhs,2);
            
            //Get x and check dimensions
            x = sizechk(prhs[1],"x",n);
            //Save x
            memcpy(objx,x,n*sizeof(double));
            
            //Create objective grad memory and get it from ASL     
            SETERRJMP();
			g = mxGetPr(plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL));
            what = "gradient";
			objgrd(0, x, g, &nerror);            
            break;
            
        case ASLCMD_CON:
            //mexPrintf("CMD: Get Constraint Val\n");
            //Check for Errors
            CHECKASL(asl);
            CHECKNRHS(nrhs,2);
            
            //Get x and check dimensions
            x = sizechk(prhs[1],"x",n);
            //Save x
            memcpy(conx,x,n*sizeof(double));
                        
            //Create constraint memory and get it from ASL  
            SETERRJMP();
			c = mxGetPr(plhs[0] = mxCreateDoubleMatrix(nc, 1, mxREAL));
            what = "constraints";
			conval(x, c, &nerror);            
            break;
            
        case ASLCMD_JAC:
            //mexPrintf("CMD: Get Jacobian Val\n");
            //Check for Errors
            CHECKASL(asl);
            CHECKNRHS(nrhs,2);
            
            //Get x and check dimensions
            x = sizechk(prhs[1],"x",n);
            //Save x
            memcpy(conx,x,n*sizeof(double));
            
            //Create constraint jac memory and get it from ASL
            SETERRJMP();            
            //Check for sparsity
            if(nrhs > 2 && *mxGetPr(prhs[2])) {
                sp = 1;
                J1 = mxGetPr(plhs[0] = mxCreateSparse(nc, n, nz, mxREAL));
            }
            else {
                sp = 0;
                J1 = mxGetPr(plhs[0] = mxCreateDoubleMatrix(nc, n, mxREAL));
            }
            
            if (nc) {
                what = "Jacobian";
                //Sparse
                if(sp) {
                    jacval(x, J1, &nerror);
                    Ir = mxGetIr(plhs[0]);
                    for(Jc = mxGetJc(plhs[0]), cs = A_colstarts, i = 0; i <= n; ++i)
                        Jc[i] = (mwIndex)cs[i];
                    cgp = Cgrad;
                    for(i = 0; i < nc; i++)
                        for(cg = *cgp++; cg; cg = cg->next)
                            Ir[cg->goff] = (mwIndex)i;  
                }
                //Dense
                else {                
                    memset(J1, 0, Jsize);
                    jacval(x, J, &nerror);
                    cgp = Cgrad;
                    for(cgpe = cgp + nc; cgp < cgpe; J1++)
                        for(cg = *cgp++; cg; cg = cg->next)
                            J1[nc*cg->varno] = J[cg->goff];
                }
            }                        
            break;
            
        case ASLCMD_JACSTR:
            //mexPrintf("CMD: Get Jacobian Structure\n");
            //Check for Errors
            CHECKASL(asl);
            CHECKNRHS(nrhs,1);
            
            //Create constraint jacstr memory and get it from ASL
            SETERRJMP();            
            
            J1 = mxGetPr(plhs[0] = mxCreateSparse(nc, n, nz, mxREAL));
            what = "Jacobian Structure)";

            for(i=0;i<nz;i++)
                J1[i] = 1.0;
            Ir = mxGetIr(plhs[0]);
            for(Jc = mxGetJc(plhs[0]), cs = A_colstarts, i = 0; i <= n; ++i)
                Jc[i] = (mwIndex)cs[i];
            cgp = Cgrad;
            for(i = 0; i < nc; i++)
                for(cg = *cgp++; cg; cg = cg->next)
                    Ir[cg->goff] = (mwIndex)i;                
            break;
            
        case ASLCMD_HES:
            //mexPrintf("CMD: Get Hessian Val\n");
            //Check for Errors
            CHECKASL(asl);
            CHECKNRHS(nrhs,4); //assume hess(x,sigma,lambda) and optionally sparse
            
            //Check dimensions & get args
            x = sizechk(prhs[1],"x",n);
            s = sizechk(prhs[2],"sigma",1);
            v = sizechk(prhs[3],"lambda",nc);
            
            //Check for sparsity
            if(nrhs > 4 && prhs[4]) {
                sp = 1;
                W = mxGetPr(plhs[0] = mxCreateSparse(n, n, nhnz, mxREAL));
            }
            else {
                sp = 0;    
                W = mxGetPr(plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL));
            }
            //Check if we need to recalculate objective / constraints
            if(!comp_x(objx,x,n))
            {
                //mexPrintf("Recalculating Objective for Hessian...\n");
                //Setup Error Catching
                SETERRJMP();
                what = "Objective for Hessian";
                //Re-evaluate Objective
                objval(0, x, &nerror);                
            }
            
            if(!comp_x(conx,x,n))
            {
                //mexPrintf("Recalculating Constraints for Hessian...\n");
                if(!c)
                    c = mxGetPr(mxCreateDoubleMatrix(nc, 1, mxREAL));
                
                //Setup Error Catching
                SETERRJMP();
                what = "Constraints for Hessian";
                //Re-evaluate Constraints
                conval(x, c, &nerror);
            }
            
            //Setup Error Catching
            SETERRJMP();
            what = "Hessian";
            
            //Sparse
            if(sp) {
                //This function returns the full (symmetric) Hessian as setup above
                sphes(H = Hsp, 1, s, v);                
                Ir = mxGetIr(plhs[0]);
                Jc = mxGetJc(plhs[0]);
                hcs = sputinfo->hcolstarts;
                hr = sputinfo->hrownos;
                for(i = 0; i <= n; i++)
                    Jc[i] = (mwIndex)hcs[i];
                He = H + hcs[n];
                while(H < He) {
                    *W++ = *H++;
                    *Ir++ = (mwIndex)*hr++;
                }	
            }
            //Dense
            else
                fullhes(W, n, 1, s, v);            
            break;
            
        case ASLCMD_HESSTR:
            //mexPrintf("CMD: Get Hessian Structure\n");
            //Check for Errors
            CHECKASL(asl);
            CHECKNRHS(nrhs,1);
            
            //Create hessianstr memory and get it from ASL
            SETERRJMP();
            W = mxGetPr(plhs[0] = mxCreateSparse(n, n, nhnz, mxREAL));
            what = "Hessian Structure";

            //sphes(H = Hsp, 0, 0, v); //i think we might have to call this...?
            Ir = mxGetIr(plhs[0]);
            Jc = mxGetJc(plhs[0]);
            hcs = sputinfo->hcolstarts;
            hr = sputinfo->hrownos;
            for(i = 0; i <= n; i++)
                Jc[i] = (mwIndex)hcs[i];
            He = H + hcs[n];
            while(H < He) {
                *W++ = 1.0;
                *H++;
                *Ir++ = (mwIndex)*hr++;
            }	                        
            break;
            
        case ASLCMD_QPHES:
            //mexPrintf("CMD: Get QP Hessian\n");
            //Check for Errors
            CHECKASL(asl);
            CHECKNRHS(nrhs,1);
            
            //Create qp hessian memory and get it from stored vals
            SETERRJMP();
            what = "QP Hessian";                         
            //Assign QP H
            if(nqpz > 0) {
                W = mxGetPr(plhs[0] = mxCreateSparse(n, n, nqpz, mxREAL));
                Ir = mxGetIr(plhs[0]);
                Jc = mxGetJc(plhs[0]);
                //Copy to Matlab Memory
                memcpy(Jc,jc,(n_var+1)*sizeof(mwIndex));
                memcpy(Ir,ir,nqpz*sizeof(mwIndex));
                memcpy(W,pr,nqpz*sizeof(double));
            }
            else
                plhs[0] = mxCreateSparse(n, n, 0, mxREAL); //empty sparse 
            break;
            
        case ASLCMD_WRITESOL:
            //Check for Errors
            CHECKASL(asl);
            CHECKNRHS(nrhs,2); //asl('writesol',msg,x)
            
            //Get Input Args
            CHECK(mxGetString(prhs[1], msg, FLEN) == 0,"error reading message!");
            x = sizechk(prhs[2],"x",n);
            
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
    else if(!strcmp(cmd,"qphes") || !strcmp(cmd,"qphess") || !strcmp(cmd,"hesqp"))
        return ASLCMD_QPHES;
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
    //Clean up any QP memory
    if(jc) {free(jc); jc = NULL;}
    if(ir) {free(ir); ir = NULL;}
    if(pr) {free(pr); pr = NULL;}
    if(objx) {free(objx); objx = NULL;}
    if(conx) {free(conx); conx = NULL;}
    
    //mexPrintf("mexExit Called\n");
    if (cur_ASL)
        ASL_free(&cur_ASL);    
}

static double* sizechk(const mxArray *mp, char *who, fint m)
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
static bool comp_x(double *xbase, double *xnew, fint n)
{
    bool same = true; fint i;
    //Loop and compare each (should be identical, no need for tolerance)
    for(i = 0; i < n; i++)
    {
        if(xbase[i] != xnew[i])
            return false;
    }
    //Default same
    return same;    
}

void printUtilityInfo()
{
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" ASL: AMPL Solver Library [v%d]\n",ASLdate_ASL);
    mexPrintf("  - Source available from: http://www.netlib.org/ampl/solvers/\n");

    mexPrintf("\n MEX Interface J.Currie 2013 (www.i2c2.aut.ac.nz)\n");
    mexPrintf("-----------------------------------------------------------\n");
}