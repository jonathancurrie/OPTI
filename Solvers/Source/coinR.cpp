/* COINR - A simple Wrapper to Use COIN-OR File Reading Routines
 * Copyright (C) 2011 Jonathan Currie (IPL)                              
 */

#include <mex.h>
#include <limits>
#include "CoinMpsIO.hpp"
#include "CoinLpIO.hpp"
#include "CoinModel.hpp"
#include "CoinMessageHandler.hpp"
#include "glpk.h"

using namespace std;

void printUtilityInfo();

//Message Handler
class DerivedHandler : public CoinMessageHandler {
public:
	virtual int print() ;
};
int DerivedHandler::print()
{
	mexPrintf(messageBuffer());
	mexPrintf("\n");
	return 0;
}
//Dummy Message Handler (for GMPL Reader - seems to crash without?)
class DummyHandler : public CoinMessageHandler {
public:
	virtual int print() ;
};
int DummyHandler::print()
{
	return 0;
}

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    //Input Args
    char *path, *type, *data = NULL;
    
    //Return Args
    double *H, *f, *A, *rl, *ru, *lb, *ub, *sostype, *sosind, *soswt;
    char *ivars;
    
    //Readers
    CoinMpsIO m;
    CoinLpIO l;
    
    //Internal Args
    mwSize buflen, dims[2] = {1,1};
    size_t i, j;
    const CoinPackedMatrix *pm;
    mwIndex *hJc, *hIr, *Jc, *Ir;
    const CoinBigIndex *jc;
    const int *ir;
    int err = -1, *ijc, *iir, printLevel = 0;
    double *elements;
    DerivedHandler *mexprinter = NULL;  
    DummyHandler *dumprinter = NULL;  
    int no_sets = 0;
    CoinSet ** sets;
    mxArray *SOSind, *SOSwt;
    int no_entries;
    const int *sos_index;
    const double *sos_weights;

    //Check Inputs
    if(nrhs < 1) {
        printUtilityInfo();        
        return;
    }   
    if(nrhs < 2) {
        mexErrMsgTxt("You must supply the file path + file type to coinR!");
        return;
    }
    if(!mxIsChar(prhs[0])) {
        mexErrMsgTxt("The file path must be a string!");
        return;
    }
    if(!mxIsChar(prhs[1])) {
        mexErrMsgTxt("The file type must be a string!");
        return;
    }
    if(nrhs > 2 && !mxIsNumeric(prhs[2])) {
        mexErrMsgTxt("Print Level must be 0 or 1!");
        return;
    }
    if(nrhs > 3 && !mxIsChar(prhs[3])) {
        mexErrMsgTxt("The data path must be a string!");
        return;
    }
    
    //File Path 
    buflen = mxGetNumberOfElements(prhs[0]) + 1; //path
    path = (char*)mxCalloc(buflen, sizeof(char)); 
    if(mxGetString(prhs[0], path, buflen) != 0)
        mexErrMsgTxt("Could not convert string data: path.");
    //File Type
    buflen = mxGetNumberOfElements(prhs[1]) + 1; //type
    type = (char*)mxCalloc(buflen, sizeof(char)); 
    if(mxGetString(prhs[1], type, buflen) != 0)
        mexErrMsgTxt("Could not convert string data: type.");
    //Print Level
    if(nrhs > 2)
        printLevel = (int)*mxGetPr(prhs[2]);  
    //Data Path
    if(nrhs > 3) {
        buflen = mxGetNumberOfElements(prhs[3]) + 1; //data
        data = (char*)mxCalloc(buflen, sizeof(char)); 
        if(mxGetString(prhs[3], type, buflen) != 0)
            mexErrMsgTxt("Could not convert string data: data.");
    }    

    //MPS or QPS or GMPL or GAMS
    if(!strcmp(type,"mps") || !strcmp(type,"qps") || !strcmp(type,"mod") || !strcmp(type,"gms")) {
        //Setup Printer
        if(printLevel) {
            mexprinter = new DerivedHandler();
            mexprinter->setLogLevel(printLevel);
            m.passInMessageHandler(mexprinter);
        }
        else if(!strcmp(type,"mod")) { //seems to crash without printer...
            dumprinter = new DummyHandler();
            dumprinter->setLogLevel(0);
            m.passInMessageHandler(dumprinter);
        }
   
        //Setup Options
        m.setInfinity(numeric_limits<double>::infinity());   
        
        //Read Based on file type
        if(!strcmp(type,"mps") || !strcmp(type,"qps"))
            err = m.readMps(path,"mps",no_sets,sets);
        else if(!strcmp(type,"mod"))
            err = m.readGMPL(path,data,true);
        else if(!strcmp(type,"gms"))
            err = m.readGms(path,"gms",true);

        if(err) {
            mexPrintf("File: %s\nError Code: %d\n",path,err);
            mexErrMsgTxt("Error Reading File! - Enable Printing to See Specific Errors");        
            return;
        }
        
        //Get Sizes
        size_t ncol = m.getNumCols();
        size_t nrow = m.getNumRows();
        if(ncol <= 0) {
            mexErrMsgTxt("Number of Columns is 0 - Error Reading File?");
            return;
        }
        //Assign ivar dims
        dims[1] = ncol;

        //Get Matrix
        pm = m.getMatrixByCol();
        size_t nelem = pm->getNumElements();
        jc = pm->getVectorStarts();
        ir = pm->getIndices();
        
        //Create Outputs    
        plhs[0] = mxCreateDoubleMatrix(ncol,1, mxREAL);     //f
        plhs[1] = mxCreateSparse(nrow,ncol,nelem,mxREAL);   //A
        Jc = mxGetJc(plhs[1]);
        Ir = mxGetIr(plhs[1]); 
        plhs[2] = mxCreateDoubleMatrix(nrow,1, mxREAL);     //rl
        plhs[3] = mxCreateDoubleMatrix(nrow,1, mxREAL);     //ru
        plhs[4] = mxCreateDoubleMatrix(ncol,1, mxREAL);     //lb
        plhs[5] = mxCreateDoubleMatrix(ncol,1, mxREAL);     //ub
        plhs[6] = mxCreateNumericArray(2,dims,mxINT8_CLASS,mxREAL); //ivars (char)
        f = mxGetPr(plhs[0]);
        A = mxGetPr(plhs[1]);
        rl = mxGetPr(plhs[2]);
        ru = mxGetPr(plhs[3]);
        lb = mxGetPr(plhs[4]);
        ub = mxGetPr(plhs[5]);
        ivars = (char*)mxGetData(plhs[6]);
        
        //Assign f
        memcpy(f,m.getObjCoefficients(),ncol*sizeof(double));
        //Assign A (Copy + Cast)
        for(i = 0; i <= ncol; i++)
            Jc[i] = (mwIndex)jc[i]; 
        for(i = 0; i < nelem; i++)
            Ir[i] = (mwIndex)ir[i];
        //Copy in A    
        memcpy(A,pm->getElements(),nelem*sizeof(double));
        //Assign Row Bounds
        if(m.getRowLower() != NULL)
            memcpy(rl,m.getRowLower(),nrow*sizeof(double));
        if(m.getRowUpper() != NULL)
            memcpy(ru,m.getRowUpper(),nrow*sizeof(double));
        //Assign Column Bounds
        if(m.getColLower() != NULL)
            memcpy(lb,m.getColLower(),ncol*sizeof(double));
        if(m.getColUpper() != NULL)
            memcpy(ub,m.getColUpper(),ncol*sizeof(double));
        //Assign integer vars
        if(m.integerColumns() != NULL)
            memcpy(ivars,m.integerColumns(),ncol*sizeof(char));    
        
        //Attempt to get QuadObj    
        err = m.readQuadraticMps(NULL, ijc, iir, elements, 0);
        if(!err) { //quadratic section read
            nelem = (size_t)ijc[ncol];
            plhs[7] = mxCreateSparse(ncol,ncol,nelem,mxREAL);   //H
            hJc = mxGetJc(plhs[7]);
            hIr = mxGetIr(plhs[7]);
            H = mxGetPr(plhs[7]);
            //Assign H (Copy + Cast)
            for(i = 0; i <= ncol; i++)
                hJc[i] = (mwIndex)ijc[i]; 
            for(i = 0; i < nelem; i++)
                hIr[i] = (mwIndex)iir[i];
            //Copy in H
            memcpy(H,elements,nelem*sizeof(double));
        }
        else //No Quad Section
            plhs[7] = mxCreateSparse(ncol,ncol,0,mxREAL); //Empty H

        //Assign problem name
        if(nlhs > 7) {
            plhs[8] = mxCreateString(m.getProblemName());
        }
        
        //Assign SOS (won't be in GMPL files)
        if(nlhs > 8) {
            plhs[9] = mxCreateDoubleMatrix(no_sets,1,mxREAL);
            plhs[10] = mxCreateCellMatrix(1,no_sets);
            plhs[11] = mxCreateCellMatrix(1,no_sets);
            if(no_sets) {
                sostype = mxGetPr(plhs[9]);
                //Assign SOS info
                for(i=0;i<no_sets;i++)
                {
                    //Assign Type
                    sostype[i] = (double)sets[i]->setType();
                    //Allocate and Assign Indices
                    no_entries = sets[i]->numberEntries();
                    sos_index = sets[i]->which();
                    SOSind = mxCreateDoubleMatrix(no_entries,1,mxREAL);
                    sosind = mxGetPr(SOSind);                    
                    for(j=0;j<no_entries;j++)
                        sosind[j] = (double)sos_index[j] + 1;
                    mxSetCell(plhs[10],i,SOSind);
                    //Allocate and Assign Weights
                    sos_weights = sets[i]->weights();
                    SOSwt = mxCreateDoubleMatrix(no_entries,1,mxREAL);
                    soswt = mxGetPr(SOSwt);                    
                    for(j=0;j<no_entries;j++)
                        soswt[j] = (double)sos_weights[j];
                    mxSetCell(plhs[11],i,SOSwt);
                }
            }
        }
        
        //Get constant objective term
        if(nlhs > 12) {
            plhs[12] = mxCreateDoubleMatrix(1,1,mxREAL);
            double objc = -m.objectiveOffset();
            memcpy(mxGetPr(plhs[12]),&objc,sizeof(double));
        }            
        
        //GMPL doesn't pick up Infinity
        if(!strcmp(type,"mod")) {
            for(i=0;i<ncol;i++){
                if(lb[i] < -1e300)
                    lb[i] = -numeric_limits<double>::infinity();
                if(ub[i] > 1e300)
                    ub[i] = numeric_limits<double>::infinity();
            }
            for(i=0;i<nrow;i++){
                if(rl[i] < -1e300)
                    rl[i] = -numeric_limits<double>::infinity();
                if(ru[i] > 1e300)
                    ru[i] = numeric_limits<double>::infinity();
            }
        }
    }
    else if(!strcmp(type,"lp"))
    {
        //Setup Printer
        if(printLevel) {
            mexprinter = new DerivedHandler();
            mexprinter->setLogLevel(printLevel);
            l.passInMessageHandler(mexprinter);
        }  
        //Setup Options
        l.setInfinity(numeric_limits<double>::infinity());
        //Read Problem
        l.readLp(path); //no error checking?

        //Get Sizes
        size_t ncol = l.getNumCols();
        size_t nrow = l.getNumRows();
        if(ncol <= 0) {
            mexErrMsgTxt("Number of Columns is 0 - Error Reading File?");
            return;
        }
        //Assign ivar dims
        dims[1] = ncol;

        //Get Matrix
        pm = l.getMatrixByCol();
        size_t nelem = pm->getNumElements();
        jc = pm->getVectorStarts();
        ir = pm->getIndices();

        //Create Outputs    
        plhs[0] = mxCreateDoubleMatrix(ncol,1, mxREAL);     //f
        plhs[1] = mxCreateSparse(nrow,ncol,nelem,mxREAL);   //A
        Jc = mxGetJc(plhs[1]);
        Ir = mxGetIr(plhs[1]); 
        plhs[2] = mxCreateDoubleMatrix(nrow,1, mxREAL);     //rl
        plhs[3] = mxCreateDoubleMatrix(nrow,1, mxREAL);     //ru
        plhs[4] = mxCreateDoubleMatrix(ncol,1, mxREAL);     //lb
        plhs[5] = mxCreateDoubleMatrix(ncol,1, mxREAL);     //ub
        plhs[6] = mxCreateNumericArray(2,dims,mxINT8_CLASS,mxREAL); //ivars (char)
        f = mxGetPr(plhs[0]);
        A = mxGetPr(plhs[1]);
        rl = mxGetPr(plhs[2]);
        ru = mxGetPr(plhs[3]);
        lb = mxGetPr(plhs[4]);
        ub = mxGetPr(plhs[5]);
        ivars = (char*)mxGetData(plhs[6]);

        //Assign f
        memcpy(f,l.getObjCoefficients(),ncol*sizeof(double));
        //Assign A (Copy + Cast)
        for(i = 0; i <= ncol; i++)
            Jc[i] = (mwIndex)jc[i]; 
        for(i = 0; i < nelem; i++)
            Ir[i] = (mwIndex)ir[i];
        //Copy in A    
        memcpy(A,pm->getElements(),nelem*sizeof(double));
        //Assign Row Bounds
        if(l.getRowLower() != NULL)
            memcpy(rl,l.getRowLower(),nrow*sizeof(double));
        if(l.getRowUpper() != NULL)
            memcpy(ru,l.getRowUpper(),nrow*sizeof(double));
        //Assign Column Bounds
        if(l.getColLower() != NULL)
            memcpy(lb,l.getColLower(),ncol*sizeof(double));
        if(l.getColUpper() != NULL)
            memcpy(ub,l.getColUpper(),ncol*sizeof(double));
        //Assign integer vars
        if(l.integerColumns() != NULL)
            memcpy(ivars,l.integerColumns(),ncol*sizeof(char));    

        //No Quad Section
        plhs[7] = mxCreateSparse(ncol,ncol,0,mxREAL); //Empty H

        //Assign problem name
        if(nlhs > 7) {
            plhs[8] = mxCreateString(l.getProblemName());
        }
        //No SOS in LP format from this reader
        if(nlhs > 8) {
            plhs[9] = mxCreateDoubleMatrix(0,1, mxREAL);  
            plhs[10] = mxCreateDoubleMatrix(0,1, mxREAL);     
            plhs[11] = mxCreateDoubleMatrix(0,1, mxREAL);  
        }
        //Get constant objective term
        if(nlhs > 12) {
            plhs[12] = mxCreateDoubleMatrix(1,1,mxREAL);
            double objc = -l.objectiveOffset();
            memcpy(mxGetPr(plhs[12]),&objc,sizeof(double));
        } 
    }
    else
        mexErrMsgTxt("Unknown Problem Type!");
    
    //Free Memory
    mxFree(path);
    mxFree(type);
}

void printUtilityInfo()
{
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" COINUTILS: COIN-OR Utilities [v%s, Built %s]\n",COINUTILS_VERSION,__DATE__);
    mexPrintf("  - Released under the Eclipse Public License: http://opensource.org/licenses/eclipse-1.0\n");
    mexPrintf("  - Source available from: https://projects.coin-or.org/CoinUtils\n\n");
    
    mexPrintf(" This binary is statically linked to the following software:\n");
    mexPrintf("  - GLPK [v%s] (GPL)\n",glp_version());
    
    mexPrintf("\n MEX Interface J.Currie 2013 (www.inverseproblem.co.nz)\n");
    mexPrintf("-----------------------------------------------------------\n");
}