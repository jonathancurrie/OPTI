/* COINW - A simple Wrapper to Use COIN-OR File Writing Routines
 * Copyright (C) 2011 Jonathan Currie (IPL)
 */

#include <mex.h>
#include <limits>
#include "CoinMpsIO.hpp"
#include "CoinLpIO.hpp"
#include "CoinModel.hpp"
#include "config_coinutils_default.h"
#include "opti_util.h"

using namespace std;
void printUtilityInfo();

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    //Input Args
    char *path, *type;
    
    //Writers
    CoinMpsIO m;
    CoinLpIO l;
    
    //Internal Args
    double *H = NULL, *f = NULL, *A = NULL, *rl = NULL, *ru = NULL, *lb = NULL, *ub = NULL;
    mwIndex *jc, *ir, *hjc, *hir;
    char *name, *ivars;
    int ncol, nrow, nelem, *Ir = NULL, *hIr = NULL, *len = NULL;
    CoinBigIndex *Jc = NULL, *hJc = NULL;
    CoinPackedMatrix *pH;
    size_t i, j;
    mwSize buflen;
    int err = -1;
    int no_sets = 0;
    CoinSet *sets = NULL;
    double *sostype, *sosind, *soswt;
    int no_entries, stype, *isosind; 
    
    //Check Inputs
    if(nrhs < 1) {
        printUtilityInfo();
        return;
    }
    if(nrhs < 3)
        mexErrMsgTxt("You must supply the problem structure + file path + file type to coinW!");
    
    if(!mxIsStruct(prhs[0]))
        mexErrMsgTxt("The problem must be a structure!");
    else {
        if(mxGetFieldNumber(prhs[0],"H") < 0)
            mexErrMsgTxt("The problem structure should contain the field 'H'");
        if(!mxIsSparse(mxGetField(prhs[0],0,"H")))
            mexErrMsgTxt("The problem objective matrix H must be sparse!'");
        if(mxGetFieldNumber(prhs[0],"f") < 0)
            mexErrMsgTxt("The problem structure should contain the field 'f'");
        if(mxIsEmpty(mxGetField(prhs[0],0,"f"))) 
            mexErrMsgTxt("The problem gradient (f) cannot be empty!");
        if(mxGetFieldNumber(prhs[0],"A") < 0)
            mexErrMsgTxt("The problem structure should contain the field 'A'");
        if(mxIsEmpty(mxGetField(prhs[0],0,"A"))) 
            mexErrMsgTxt("The problem constraint matrix (A) cannot be empty!");
        if(!mxIsSparse(mxGetField(prhs[0],0,"A")))
            mexErrMsgTxt("The problem constraint matrix A must be sparse!'");
        if(mxGetFieldNumber(prhs[0],"rl") < 0)
            mexErrMsgTxt("The problem structure should contain the field 'rl'");
        if(mxGetFieldNumber(prhs[0],"ru") < 0)
            mexErrMsgTxt("The problem structure should contain the field 'ru'");
        if(mxGetFieldNumber(prhs[0],"lb") < 0)
            mexErrMsgTxt("The problem structure should contain the field 'lb'");
        if(mxGetFieldNumber(prhs[0],"ub") < 0)
            mexErrMsgTxt("The problem structure should contain the field 'ub'");
        if(mxGetFieldNumber(prhs[0],"int") < 0)
            mexErrMsgTxt("The problem structure should contain the field 'int'");
        if(mxGetFieldNumber(prhs[0],"Name") < 0)
            mexErrMsgTxt("The problem structure should contain the field 'Name'");
        if(mxGetFieldNumber(prhs[0],"sos_type") < 0)
            mexErrMsgTxt("The problem structure should contain the field 'sos_type'");
        if(mxGetFieldNumber(prhs[0],"sos_index") < 0)
            mexErrMsgTxt("The problem structure should contain the field 'sos_index'");
        if(mxGetFieldNumber(prhs[0],"sos_weight") < 0)
            mexErrMsgTxt("The problem structure should contain the field 'sos_weight'");
    }  
    if(!mxIsChar(prhs[1]))
        mexErrMsgTxt("The file path must be a string!");
    if(!mxIsChar(prhs[2]))
        mexErrMsgTxt("The file type must be a string!");
    if(nrhs > 3 && !mxIsNumeric(prhs[3]))
        mexErrMsgTxt("Print Level must be 0 or 1!");
    
    //Get Problem Elements
    if(!mxIsEmpty(mxGetField(prhs[0],0,"H"))) {
        H = mxGetPr(mxGetField(prhs[0],0,"H"));
        hjc = mxGetJc(mxGetField(prhs[0],0,"H"));
        hir = mxGetIr(mxGetField(prhs[0],0,"H"));
    }
    f = mxGetPr(mxGetField(prhs[0],0,"f"));
    A = mxGetPr(mxGetField(prhs[0],0,"A"));
    jc = mxGetJc(mxGetField(prhs[0],0,"A"));
    ir = mxGetIr(mxGetField(prhs[0],0,"A"));
    rl = mxGetPr(mxGetField(prhs[0],0,"rl"));
    ru = mxGetPr(mxGetField(prhs[0],0,"ru"));
    if(!mxIsEmpty(mxGetField(prhs[0],0,"lb")))
        lb = mxGetPr(mxGetField(prhs[0],0,"lb"));
    if(!mxIsEmpty(mxGetField(prhs[0],0,"ub")))
        ub = mxGetPr(mxGetField(prhs[0],0,"ub"));
    ivars = (char*)mxGetData(mxGetField(prhs[0],0,"int"));
    
    //Get Sizes
    ncol = (int)mxGetNumberOfElements(mxGetField(prhs[0],0,"f"));
    nrow = (int)mxGetM(mxGetField(prhs[0],0,"A"));
    
    //File Path 
    buflen = mxGetNumberOfElements(prhs[1]) + 1; //path
    path = (char*)mxCalloc(buflen, sizeof(char)); 
    if(mxGetString(prhs[1], path, buflen) != 0)
        mexErrMsgTxt("Could not convert string data: path.");
    //File Type
    buflen = mxGetNumberOfElements(prhs[2]) + 1; //type
    type = (char*)mxCalloc(buflen, sizeof(char)); 
    if(mxGetString(prhs[2], type, buflen) != 0)
        mexErrMsgTxt("Could not convert string data: type.");

    //Build A (copy + cast)
    nelem = (int)jc[ncol];
    Jc = (CoinBigIndex*)mxCalloc(ncol+1,sizeof(CoinBigIndex));
    Ir = (int*)mxCalloc(nelem,sizeof(int));
    for(i=0;i<=ncol;i++)
        Jc[i] = (CoinBigIndex)jc[i];
    for(i=0;i<nelem;i++)
        Ir[i] = (int)ir[i];
    //Assign to packed matrix  
    CoinPackedMatrix pm(true,nrow,ncol,nelem,A,Ir,Jc,len);
    
    
    if(!strcmp(type,"mps") || !strcmp(type,"qps")) {  
        //Build H (copy + cast)
        if(H) {
            nelem = (int)jc[ncol];
            hJc = (CoinBigIndex*)mxCalloc(ncol+1,sizeof(CoinBigIndex));
            hIr = (int*)mxCalloc(nelem,sizeof(int));
            for(i=0;i<=ncol;i++)
                hJc[i] = (CoinBigIndex)hjc[i];
            for(i=0;i<nelem;i++)
                hIr[i] = (int)hir[i];
            //Assign to packed matrix  
            pH = new CoinPackedMatrix(true,nrow,ncol,nelem,H,hIr,hJc,len);
        }

        //Set Problem Data
        m.setMpsData(pm,numeric_limits<double>::infinity(),lb,ub,f,ivars,rl,ru,NULL,NULL);
        //Set Name
        buflen = mxGetNumberOfElements(mxGetField(prhs[0],0,"Name")) + 1; //type
        name = (char*)mxCalloc(buflen, sizeof(char)); 
        if(mxGetString(mxGetField(prhs[0],0,"Name"), name, buflen) != 0)
            mexErrMsgTxt("Could not convert string data: name.");
        m.setProblemName(name) ;               
        
        //Set SOS
        no_sets = (int)mxGetNumberOfElements(mxGetField(prhs[0],0,"sos_type"));
        if(no_sets > 0) {
            //Error Checking
            if(!mxIsCell(mxGetField(prhs[0],0,"sos_index")) || mxIsEmpty(mxGetField(prhs[0],0,"sos_index")))
                mexErrMsgTxt("sos_index must be a cell array, and not empty!");
            if(!mxIsCell(mxGetField(prhs[0],0,"sos_weight")) || mxIsEmpty(mxGetField(prhs[0],0,"sos_weight")))
                mexErrMsgTxt("sos_weight must be a cell array, and not empty!");
            if(mxGetNumberOfElements(mxGetField(prhs[0],0,"sos_index")) != no_sets)
                mexErrMsgTxt("sos_index cell array is not the same length as sos_type!");
            if(mxGetNumberOfElements(mxGetField(prhs[0],0,"sos_weight")) != no_sets)
                mexErrMsgTxt("sos_weight cell array is not the same length as sos_type!");
            
            //Collect Types
            sostype = mxGetPr(mxGetField(prhs[0],0,"sos_type"));
            
            //Allocate Set Memory
            sets = new CoinSet[no_sets];
            //Copy in SOS data, creating CoinSosSets as we go
            for(i=0;i<no_sets;i++) {
                stype = (int)sostype[i];
                no_entries = (int)mxGetNumberOfElements(mxGetCell(mxGetField(prhs[0],0,"sos_index"),i));
                //Create sosind memory and copy in
                isosind = new int[no_entries];
                sosind = mxGetPr(mxGetCell(mxGetField(prhs[0],0,"sos_index"),i));
                for(j=0;j<no_entries;j++)
                    isosind[j] = (int)sosind[j]-1;                
                //Get soswt
                soswt = mxGetPr(mxGetCell(mxGetField(prhs[0],0,"sos_weight"),i));
                //Create Set object
                sets[i] = CoinSosSet(no_entries,isosind,soswt,stype);
                
                delete isosind; //free memory as we go round, CoinSet copies internally
            }
        }
        
        //Set constant objective term
        if(mxGetField(prhs[0],0,"objbias") && !mxIsEmpty(mxGetField(prhs[0],0,"objbias")))
             m.setObjectiveOffset(-*mxGetPr(mxGetField(prhs[0],0,"objbias")));                  
        
        //Write Problem
        if(H)
            err = m.writeMps(path, 0, 1, 2, pH, no_sets, sets); //as below but with quad
        else
            err = m.writeMps(path, 0, 1, 2, NULL, no_sets, sets); //no compression, enhanced precision, 2 vals across, no quad
    }
    else if(!strcmp(type,"lp")) {
        if(H)
            mexErrMsgTxt("You cannot write a QP to an LP file!");

        no_sets = (int)mxGetNumberOfElements(mxGetField(prhs[0],0,"sos_type"));
        if(no_sets > 0) 
            mexWarnMsgTxt("Currently SOS cannot be written to a LP file");
        
        //Set Problem Data
        l.setLpDataWithoutRowAndColNames(pm,lb,ub,f,ivars,rl,ru);
        //Set Name
        buflen = mxGetNumberOfElements(mxGetField(prhs[0],0,"Name")) + 1; //type
        name = (char*)mxCalloc(buflen, sizeof(char)); 
        if(mxGetString(mxGetField(prhs[0],0,"Name"), name, buflen) != 0)
            mexErrMsgTxt("Could not convert string data: name.");
        l.setProblemName(name) ;
        
        //Set constant objective term
        if(mxGetField(prhs[0],0,"objbias") && !mxIsEmpty(mxGetField(prhs[0],0,"objbias")))
             l.setObjectiveOffset(-*mxGetPr(mxGetField(prhs[0],0,"objbias"))); 

        //Write Problem
        err = l.writeLp(path); //use Row names
    }
    else
        mexErrMsgTxt("Unknown problem type!");
    
    //Free Memory
    mxFree(Jc);
    mxFree(Ir);
    mxFree(path);
    mxFree(type);
    mxFree(name);
}

void printUtilityInfo()
{
    char vbuf[6]; getVSVer(vbuf); 
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" COINUTILS: COIN-OR Utilities [v%s, Built %s, VS%s]\n",COINUTILS_VERSION,__DATE__,vbuf);
    mexPrintf("  - Released under the Eclipse Public License: http://opensource.org/licenses/eclipse-1.0\n");
    mexPrintf("  - Source available from: https://projects.coin-or.org/CoinUtils\n");

    mexPrintf("\n MEX Interface J.Currie 2013 (www.inverseproblem.co.nz)\n");
    mexPrintf("-----------------------------------------------------------\n");
}