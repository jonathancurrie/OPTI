#ifndef INCLUDE_LBFGSB
#define INCLUDE_LBFGSB

//Function handle structure
#define FLEN 128 /* max length of user function name */
#define MAXRHS 2 /* max nrhs for user function */
typedef struct {
     char f[FLEN];
     mxArray *plhs[1];
     mxArray *prhs[MAXRHS];
     int xrhs, nrhs;
} user_function_data;

//Iteration callback structure
typedef struct {
    char f[FLEN];
    mxArray *plhs[1];
    mxArray *prhs[3];
    bool enabled;
} iter_fun_data;

//Macros
#define CHECK(cond, msg) if (!(cond)) { mexErrMsgTxt(msg); }

//Ctrl-C Detection (Undocumented - Found in gurobi_mex.c!)
#ifdef __cplusplus
    extern "C" bool utIsInterruptPending();
    extern "C" void utSetInterruptPending(bool);
#else
    extern bool utIsInterruptPending();
    extern void utSetInterruptPending(bool);
#endif

#endif