#include "mex.h"

#include "lp_lib.h"

#if (defined WIN32 || defined _WIN32 || defined MSDOS || defined DOS) && !defined LPSOLVEAPIFROMLIB
# define LPSOLVEAPIEXPLICIT

# include "lp_explicit.h"
#else
# undef LPSOLVEAPIFROMLIB

#if 1
/* static linking */
# define init_lpsolve_lib() TRUE
# define exit_lpsolve_lib()
# define putlogfunc   put_logfunc
# define putabortfunc put_abortfunc

#else
/* dynamic linking */
#if defined driverVERSION

# define LPSOLVEAPIFROMLIB

# include "lp_explicit.h"

static hlpsolve hlpsolve_ = NULL;

#endif

# define init_lpsolve_lib() ((hlpsolve_ != NULL) || (((hlpsolve_ = open_lpsolve_lib(NULL)) != NULL) && init_lpsolve(hlpsolve_)))
# define exit_lpsolve_lib() { if (hlpsolve_ != NULL) close_lpsolve_lib(hlpsolve_); hlpsolve_ = NULL; }
# define putlogfunc put_logfunc
# define putabortfunc put_abortfunc
#endif

#endif

#define quotechar "'"
#define ErrMsgTxt(lpsolvecaller, str) mexErrMsgTxt(str)
#define drivername mxlpsolve
#define strdrivername "lp_solve" //modified
#define caller "MATLAB"
#define IsNumeric mxIsNumeric
#define IsComplex mxIsComplex
#define IsSparse mxIsSparse
#define GetM(lpsolvecaller, mat) mxGetM(mat)
#define GetN(lpsolvecaller, mat) mxGetN(mat)
#define GetpMatrix(lpsolvecaller, element) (lpsolvecaller)->prhs[element]
#define GetPr mxGetPr
#define matCalloc mxCalloc
#define matRealloc mxRealloc
#define matFree mxFree
#define GetCellString(lpsolvecaller, ppm, element, buf, size) GetString(lpsolvecaller, ppm, element, buf, size, TRUE)
#define Printf mexPrintf

#define MatrixEl mxArray *
#define pMatrix MatrixEl *
#define rMatrix MatrixEl
#define strArray pMatrix

#define callerPrototype(callername) void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])

#define publicargs(lpsolve) \
        (lpsolve)->lpsolvecaller.plhs = (pMatrix) plhs; \
        (lpsolve)->lpsolvecaller.prhs = (pMatrix) prhs; \
        (lpsolve)->lpsolvecaller.nlhs = nlhs; \
        (lpsolve)->lpsolvecaller.nrhs = nrhs;

#define registerExitFcn(lpsolve) if (mexAtExit(ExitFcn)) ErrMsgTxt(&lpsolve->lpsolvecaller, "Failed to register exit function.\n")

#define ExitcallerPrototype(lpsolve) return

#define BEGIN_INTERRUPT_IMMEDIATELY_IN_FOREIGN_CODE
#define END_INTERRUPT_IMMEDIATELY_IN_FOREIGN_CODE


typedef struct
{
        jmp_buf exit_mark;
        int nlhs;
        int nrhs;
        pMatrix plhs;
        pMatrix prhs;
} structlpsolvecaller;

#define Double double
#define Long Double

#define CreateLongMatrix CreateDoubleMatrix
#define SetDoubleMatrix(lpsolvecaller, mat, m, n, element, freemat)
#define SetLongMatrix SetDoubleMatrix
