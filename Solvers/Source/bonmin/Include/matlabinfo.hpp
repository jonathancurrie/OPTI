// Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         September 25, 2008

// Modified J.Currie September 2011 to suit BONMIN

#ifndef INCLUDE_MATLABINFO
#define INCLUDE_MATLABINFO

#include "mex.h"

// Class MatlabInfo.
// -----------------------------------------------------------------
// An object of this class stores all the information we will pass 
// back to MATLAB upon termination of IPOPT.
class MatlabInfo {
public:

  // Create a new info object and store the information in a MATLAB
  // array. The input pointer will point to the the newly created
  // MATLAB array. Since the user has an opportunity to modify the
  // MATLAB array pointed to by "ptr", we do not destroy the array
  // when the MatlabInfo object is destroyed. It is up to the user to
  // do that.
  explicit MatlabInfo (mxArray*& ptr);

  // The destructor.
  ~MatlabInfo() { };

  // Access and modify the exit status and solution statistics.
  int       getExitStatus () const;
  void      setExitStatus (int status);
  void 	    setNumNodes (int nodes);
  void      setIterationCount (int iter);
  void 	    setBestObj (double obj);
  void      setCpuTime (double cpu);

protected:
  mxArray* ptr;  // All the information is stored in a MATLAB array.
};

#endif
