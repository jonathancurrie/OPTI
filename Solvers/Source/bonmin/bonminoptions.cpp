// Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         September 15, 2008

// Modified J.Currie September 2011 + Jan 2013 to suit BONMIN

#include "bonminoptions.hpp"
#include "matlabexception.hpp"

 using Ipopt::SmartPtr;
 using Ipopt::IsValid;
 using Ipopt::RegisteredOption;
 using Ipopt::Snprintf;


// Function definitions for class BonminOptions.
// -----------------------------------------------------------------
BonminOptions::BonminOptions (BonminSetup& app, const mxArray* ptr) 
  : app(app) {

    // Check to make sure the MATLAB array is a structure array. 
    if (!mxIsStruct(ptr))
    throw MatlabException("The OPTIONS input must be a structure array! Type HELP STRUCT in the MATLAB console for more information");

    //Check for IPOPT options, process if we found them
    mxArray *iOpt = mxGetField(ptr,0,"ipopt");
    if(iOpt) {
        // Each field in the structure array should correspond to an option in IPOPT. Repeat for each field.
        int n = mxGetNumberOfFields(iOpt);
        for (int i = 0; i < n; i++) {
            const char* label = mxGetFieldNameByNumber(iOpt,i);
            mxArray*    p     = mxGetFieldByNumber(iOpt,0,i);
            setOption(label,p);            
        }
    }
    //Check for BONMIN options, process if we found them
    mxArray *bOpt = mxGetField(ptr,0,"bonmin");
    if(bOpt) {
        // Each field in the structure array should correspond to an option in BONMIN. Repeat for each field.
        int n = mxGetNumberOfFields(bOpt);
        for (int i = 0; i < n; i++) {
            const char* label = mxGetFieldNameByNumber(bOpt,i);
            mxArray*    p     = mxGetFieldByNumber(bOpt,0,i);  
            setOption(label,p);
        }
    }
}

bool BonminOptions::useQuasiNewton() const {
  bool        b;  // The return value.
  std::string value;

  app.options()->GetStringValue("hessian_approximation",value,"");
  b = !value.compare("limited-memory");
  return b;
}

bool BonminOptions::useDerivChecker() const {
  bool        b;  // The return value.
  std::string value;
  app.options()->GetStringValue("derivative_test",value,"");  
  b = value.compare("none");
  return b;
}

bool BonminOptions::userScaling() const {
  bool        b;  // The return value.
  std::string value;
  app.options()->GetStringValue("nlp_scaling_method",value,"");  
  b = !value.compare("user-scaling");
  return b;
}

int BonminOptions::printLevel() const {
  int value;  // The return value.
  app.options()->GetIntegerValue("print_level",value,"");
  return value;
}

bool BonminOptions::usingMA57() const {
  bool b;
  std::string value;

  app.options()->GetStringValue("linear_solver",value,"");  
  b = value.compare("ma57");
  return b;
}

double BonminOptions::getPosInfty() const {
  //double value;  // The return value.
  //app.options()->GetNumericValue("nlp_upper_bound_inf",value,"");
  return DBL_MAX;
}

double BonminOptions::getNegInfty() const {
  //double value;  // The return value.
  //app.options()->GetNumericValue("nlp_lower_bound_inf",value,"");
  return -DBL_MAX;
}

void BonminOptions::setOption (const char* label, const mxArray* ptr) {

  // Check to make sure we have a valid option.
  SmartPtr<const RegisteredOption> option = app.roptions()->GetOption(label);
  if (!IsValid(option)) {
    char buf[256];
    Snprintf(buf, 255, "You have specified a nonexistent BONMIN/IPOPT option (\"%s\")", label);
    throw MatlabException(buf);
  }

  Ipopt::RegisteredOptionType type = option->Type();
  if (type == Ipopt::OT_String)
    setStringOption(label,ptr);
  else if (type == Ipopt::OT_Integer)
    setIntegerOption(label,ptr);
  else
    setNumberOption(label,ptr);
}

void BonminOptions::printOption(const char *label)
{
    // Check to make sure we have a valid option.
    SmartPtr<const RegisteredOption> option = app.roptions()->GetOption(label);
    if (!IsValid(option)) {
        char buf[256];
        Snprintf(buf, 255, "You have specified a nonexistent BONMIN/IPOPT option (\"%s\")", label);
        throw MatlabException(buf);
    }
    Ipopt::RegisteredOptionType type = option->Type();
    if (type == Ipopt::OT_Integer) {
        int val;
        app.options()->GetIntegerValue(label,val,"");
        mexPrintf("Integer Option %s = %d\n",label,val);
    }
    else if (type == Ipopt::OT_String) {
        std::string val;
        app.options()->GetStringValue(label,val,"");  
        mexPrintf("String Option %s = %s\n",label,val);
    }
    else {
        double val;
        app.options()->GetNumericValue(label,val,"");  
        mexPrintf("Real Option %s = %f\n",label,val);
    } 
}

void BonminOptions::setStringOption (const char* label, const mxArray* ptr) {

  // Check whether the option value is a string.
  if (!mxIsChar(ptr)) {
    char buf[256];
    Snprintf(buf, 255, "BONMIN/IPOPT option value for option \"%s\" should be a string", label);
    throw MatlabException(buf);
  }

  // Get the option value.
  char* value = mxArrayToString(ptr);

  // Set the option.
  bool success = app.options()->SetStringValue(label,value);
  if (!success) {
    char buf[256];
    Snprintf(buf, 255, "Invalid value for BONMIN/IPOPT option \"%s\"", label);
    throw MatlabException(buf);
  }

  // Free the dynamically allocated memory.
  mxFree(value);
}

void BonminOptions::setIntegerOption (const char* label, const mxArray* ptr) {
  
  // Check whether the option value is a number.
  if (!mxIsDouble(ptr)) {
    char buf[256];
    Snprintf(buf, 255, "BONMIN/IPOPT option value for option \"%s\" should be an integer", label);
    throw MatlabException(buf);
  }
  
  // Set either the integer option.
  double value   = mxGetScalar(ptr);
  bool   success = app.options()->SetIntegerValue(label,(int) value);
  if (!success) {
    char buf[256];
    Snprintf(buf, 255, "Invalid value for integer BONMIN/IPOPT option \"%s\"", label);
    throw MatlabException(buf);
  }
}

void BonminOptions::setNumberOption (const char* label, const mxArray* ptr) {
  
  // Check whether the option value is a number.
  if (!mxIsDouble(ptr)) {
    char buf[256];
    Snprintf(buf, 255, "BONMIN/IPOPT option value for option \"%s\" should be a number", label);
    throw MatlabException(buf);
  }
  
  // Set either the numeric option.
  double value   = mxGetScalar(ptr);
  bool   success = app.options()->SetNumericValue(label,value);
  if (!success) {
    char buf[256];
    Snprintf(buf, 255, "Invalid value for numeric BONMIN/IPOPT option \"%s\"", label);
    throw MatlabException(buf);
  }
}
