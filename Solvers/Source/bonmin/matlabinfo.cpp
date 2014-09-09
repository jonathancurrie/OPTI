// Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         September 25, 2008

// Modified J.Currie September 2011 to suit BONMIN

#include "matlabinfo.hpp"
#include "iterate.hpp"

// Function definitions for class MatlabInfo.
// ------------------------------------------------------------------
MatlabInfo::MatlabInfo (mxArray*& ptr) 
  : ptr(0) {

  // Create the structure array.
  const char* fieldnames[5];
  const char* exitstatusfield = "status";
  const char* iterfield       = "iter";
  const char* nodesfield      = "nodes";
  const char* bestfield       = "bestobj";
  const char* cpu             = "cpu";
  fieldnames[0] = exitstatusfield;
  fieldnames[1] = iterfield;
  fieldnames[2] = nodesfield;
  fieldnames[3] = bestfield;
  fieldnames[4] = cpu;
  this->ptr = ptr = mxCreateStructMatrix(1,1,5,fieldnames);

  // Initialize some fields.
  mxSetField(ptr,0,"status",mxCreateDoubleScalar(0));
  mxSetField(ptr,0,"iter",mxCreateDoubleScalar(0));
  mxSetField(ptr,0,"nodes",mxCreateDoubleScalar(0));
  mxSetField(ptr,0,"bestobj",mxCreateDoubleScalar(0));
  mxSetField(ptr,0,"cpu",mxCreateDoubleScalar(0));
}

int MatlabInfo::getExitStatus() const {
  const mxArray* p = mxGetField(ptr,0,"status");
  return (int) *mxGetPr(p);  
}

void MatlabInfo::setExitStatus (int status) {
  mxArray* p = mxGetField(ptr,0,"status");
  *mxGetPr(p) = (double) status;
}

void MatlabInfo::setNumNodes (int nodes) {
  mxArray* p = mxGetField(ptr,0,"nodes");
  *mxGetPr(p) = (double) nodes;
}

void MatlabInfo::setIterationCount (int iter) {
  mxArray* p = mxGetField(ptr,0,"iter");
  *mxGetPr(p) = (double) iter;
}

void MatlabInfo::setBestObj (double obj) {
  mxArray* p = mxGetField(ptr,0,"bestobj");
  *mxGetPr(p) = obj;
}

void MatlabInfo::setCpuTime (double cpu) {
  mxArray* p = mxGetField(ptr,0,"cpu");
  *mxGetPr(p) = cpu;
}
