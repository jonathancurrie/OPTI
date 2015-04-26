/*
    Author:        Peter Notebaert
    Contact:       lpsolve@peno.be
    License terms: LGPL.

    Release notes:
      interface to the lp_solve 5.5 toolkit.

      This file and the MATLAB interface code is originally developed by
      Jeffrey C. Kantor (jeff@control.cheg.nd.edu, cchen@darwin.helion.nd.edu)
      for lp_solve version 2.3
      The original name was lpmex.c


      It is completely revised and redesigned by Peter Notebaert for lp_solve version 5.5

      lpsolve.c needs to be linked with hash.c and the lp_solve library.

      The hash function is needed to have a quick access from the provided functionname to the
      implementation of the function.

    Change history:
     v5.1.0.0
      - First implementation
     v5.5.0.1  12 nov 2005
      - set_outputfile(lp0, ""); added to disable messages on console
      - For routines that return an integer,
        CreateDoubleMatrix and SetDoubleMatrix replaced by CreateLongMatrix and SetLongMatrix
     v5.5.0.2  19 nov 2005
      - For routines that need an lp handle, this handle can now also be the model name.
      - New call get_handle to get handle from model name.
     v5.5.0.4  25 aug 2006
      - New routine guess_basis added.
     v5.5.0.5  31 okt 2006
      - Fixed a bug in hash routines which could result in a severe error with
        model names instead of handles
     v5.5.0.6  14 okt 2008
      - Revised for PHP
      - Removed static variables for multithreading.
     v5.5.0.7  4 may 2009
      - Sometimes a memory leak could occur. Fixed.
      - Added the possibility to call a check routine in the abort callback routine.
      - Added the possibility to enable fortify (just add -DFORTIFY at compile command line)
     v5.5.0.8  26 jul 2009
      - Constants can now also be provided as string. For example 'LE'
      - Returned constants can be returned as string by calling return_constants with a non-zero parameter
        When return_constants is called without a parameter, it only returns the current setting
*/

/* Modified J.Currie September 2011 */

#include <stdio.h>
#include <signal.h>
#include <string.h>
#include "opti_util.h"

#include "lpsolvecaller.h"

#define NEWLINE "\n"

/* Declare a global lp record */

#define LPSTEP 100

void printSolverInfo();

static lprec   **lp;
static hashtable *cmdhash, *constanthash, *handlehash;
static int     lp_last;
static short   initialized = FALSE;
static short   interrupted = FALSE;
static char    return_constants = FALSE; /* return constants as string N/Y */

#define bufsz 200

enum consttypes { /* must be a power of 2 */
 consttype_constrainttype = 1,
 consttype_antidegen = 2,
 consttype_branch = 4,
 consttype_crash = 8,
 consttype_verbose = 16,
 consttype_solve = 32,
 consttype_improve = 64,
 consttype_msg = 128,
 consttype_node = 256,
 consttype_presolve = 512,
 consttype_pricer = 1024,
 consttype_price = 2048,
 consttype_scale = 4096,
 consttype_simplex = 8192,
 consttype_eps = 16384,
 consttype_MPS = 32768,
};

static struct {
        char *svalue;
        int value;
        unsigned int consttype;   /* one of the consttypes enums. It can be an or of multiple of them */
        unsigned int mask;        /* if 0, then value is a bit mask. Else value is a value when masked with this value */
        unsigned char usereturn;  /* use this constant as return value. If not then it is only used as input constant */
} constants[] = {
  { "<", LE, consttype_constrainttype, ~0, FALSE },
  { "<=", LE, consttype_constrainttype, ~0, FALSE },
  { ">", GE, consttype_constrainttype, ~0, FALSE },
  { ">=", GE, consttype_constrainttype, ~0, FALSE },
  { "=", EQ, consttype_constrainttype, ~0, FALSE },
  { "FR", FR, consttype_constrainttype, ~0, TRUE },
  { "LE", LE, consttype_constrainttype, ~0, TRUE },
  { "GE", GE, consttype_constrainttype, ~0, TRUE },
  { "EQ", EQ, consttype_constrainttype, ~0, TRUE },

  { "NEUTRAL", NEUTRAL, consttype_verbose, 0x07, TRUE },
  { "CRITICAL", CRITICAL, consttype_verbose, 0x07, TRUE },
  { "SEVERE", SEVERE, consttype_verbose, 0x07, TRUE },
  { "IMPORTANT", IMPORTANT, consttype_verbose, 0x07, TRUE },
  { "NORMAL", NORMAL, consttype_verbose, 0x07, TRUE },
  { "DETAILED", DETAILED, consttype_verbose, 0x07, TRUE },
  { "FULL", FULL, consttype_verbose, 0x07, TRUE },

  { "MPS_FREE", MPS_FREE, consttype_MPS, 0, TRUE },
  { "MPS_IBM", MPS_IBM, consttype_MPS, 0, TRUE },
  { "MPS_NEGOBJCONST", MPS_NEGOBJCONST, consttype_MPS, 0, TRUE },

  { "NOMEMORY", NOMEMORY, consttype_solve, ~0, TRUE },
  { "OPTIMAL", OPTIMAL, consttype_solve, ~0, TRUE },
  { "SUBOPTIMAL", SUBOPTIMAL, consttype_solve, ~0, TRUE },
  { "INFEASIBLE", INFEASIBLE, consttype_solve, ~0, TRUE },
  { "UNBOUNDED", UNBOUNDED, consttype_solve, ~0, TRUE },
  { "DEGENERATE", DEGENERATE, consttype_solve, ~0, TRUE },
  { "NUMFAILURE", NUMFAILURE, consttype_solve, ~0, TRUE },
  { "USERABORT", USERABORT, consttype_solve, ~0, TRUE },
  { "TIMEOUT", TIMEOUT, consttype_solve, ~0, TRUE },
  { "PRESOLVED", PRESOLVED, consttype_solve, ~0, TRUE },
  { "PROCFAIL", PROCFAIL, consttype_solve, ~0, TRUE },
  { "PROCBREAK", PROCBREAK, consttype_solve, ~0, TRUE },
  { "FEASFOUND", FEASFOUND, consttype_solve, ~0, TRUE },
  { "NOFEASFOUND", NOFEASFOUND, consttype_solve, ~0, TRUE },

  { "ANTIDEGEN_NONE", ANTIDEGEN_NONE, consttype_antidegen, ~0, TRUE },
  { "ANTIDEGEN_FIXEDVARS", ANTIDEGEN_FIXEDVARS, consttype_antidegen, 0, TRUE },
  { "ANTIDEGEN_COLUMNCHECK", ANTIDEGEN_COLUMNCHECK, consttype_antidegen, 0, TRUE },
  { "ANTIDEGEN_STALLING", ANTIDEGEN_STALLING, consttype_antidegen, 0, TRUE },
  { "ANTIDEGEN_NUMFAILURE", ANTIDEGEN_NUMFAILURE, consttype_antidegen, 0, TRUE },
  { "ANTIDEGEN_LOSTFEAS", ANTIDEGEN_LOSTFEAS, consttype_antidegen, 0, TRUE },
  { "ANTIDEGEN_INFEASIBLE", ANTIDEGEN_INFEASIBLE, consttype_antidegen, 0, TRUE },
  { "ANTIDEGEN_DYNAMIC", ANTIDEGEN_DYNAMIC, consttype_antidegen, 0, TRUE },
  { "ANTIDEGEN_DURINGBB", ANTIDEGEN_DURINGBB, consttype_antidegen, 0, TRUE },
  { "ANTIDEGEN_RHSPERTURB", ANTIDEGEN_RHSPERTURB, consttype_antidegen, 0, TRUE },
  { "ANTIDEGEN_BOUNDFLIP", ANTIDEGEN_BOUNDFLIP, consttype_antidegen, 0, TRUE },

  { "BRANCH_CEILING", BRANCH_CEILING, consttype_branch, ~0, TRUE },
  { "BRANCH_FLOOR", BRANCH_FLOOR, consttype_branch, ~0, TRUE },
  { "BRANCH_AUTOMATIC", BRANCH_AUTOMATIC, consttype_branch, ~0, TRUE },
  { "BRANCH_DEFAULT", BRANCH_DEFAULT, consttype_branch, ~0, TRUE },

  { "CRASH_NONE", CRASH_NONE, consttype_crash, ~0, TRUE },
  { "CRASH_MOSTFEASIBLE", CRASH_MOSTFEASIBLE, consttype_crash, ~0, TRUE },
  { "CRASH_LEASTDEGENERATE", CRASH_LEASTDEGENERATE, consttype_crash, ~0, TRUE },

  { "IMPROVE_NONE", IMPROVE_NONE, consttype_improve, ~0, TRUE },
  { "IMPROVE_SOLUTION", IMPROVE_SOLUTION, consttype_improve, 0, TRUE },
  { "IMPROVE_DUALFEAS", IMPROVE_DUALFEAS, consttype_improve, 0, TRUE },
  { "IMPROVE_THETAGAP", IMPROVE_THETAGAP, consttype_improve, 0, TRUE },
  { "IMPROVE_BBSIMPLEX", IMPROVE_BBSIMPLEX, consttype_improve, 0, TRUE },

  { "MSG_PRESOLVE", MSG_PRESOLVE, consttype_msg, 0, TRUE },
  { "MSG_LPFEASIBLE", MSG_LPFEASIBLE, consttype_msg, 0, TRUE },
  { "MSG_LPOPTIMAL", MSG_LPOPTIMAL, consttype_msg, 0, TRUE },
  { "MSG_MILPFEASIBLE", MSG_MILPFEASIBLE, consttype_msg, 0, TRUE },
  { "MSG_MILPEQUAL", MSG_MILPEQUAL, consttype_msg, 0, TRUE },
  { "MSG_MILPBETTER", MSG_MILPBETTER, consttype_msg, 0, TRUE },

  { "NODE_FIRSTSELECT", NODE_FIRSTSELECT, consttype_node, 0x07, TRUE },
  { "NODE_GAPSELECT", NODE_GAPSELECT, consttype_node, 0x07, TRUE },
  { "NODE_RANGESELECT", NODE_RANGESELECT, consttype_node, 0x07, TRUE },
  { "NODE_FRACTIONSELECT", NODE_FRACTIONSELECT, consttype_node, 0x07, TRUE },
  { "NODE_PSEUDOCOSTSELECT", NODE_PSEUDOCOSTSELECT, consttype_node, 0x07, TRUE },
  { "NODE_PSEUDONONINTSELECT", NODE_PSEUDONONINTSELECT, consttype_node, 0x07, TRUE },
  { "NODE_PSEUDORATIOSELECT", NODE_PSEUDORATIOSELECT, consttype_node, 0x07, TRUE },
  { "NODE_USERSELECT", NODE_USERSELECT, consttype_node, 0x07, TRUE },
  { "NODE_WEIGHTREVERSEMODE", NODE_WEIGHTREVERSEMODE, consttype_node, 0, TRUE },
  { "NODE_BRANCHREVERSEMODE", NODE_BRANCHREVERSEMODE, consttype_node, 0, TRUE },
  { "NODE_GREEDYMODE", NODE_GREEDYMODE, consttype_node, 0, TRUE },
  { "NODE_PSEUDOCOSTMODE", NODE_PSEUDOCOSTMODE, consttype_node, 0, TRUE },
  { "NODE_DEPTHFIRSTMODE", NODE_DEPTHFIRSTMODE, consttype_node, 0, TRUE },
  { "NODE_RANDOMIZEMODE", NODE_RANDOMIZEMODE, consttype_node, 0, TRUE },
  { "NODE_GUBMODE", NODE_GUBMODE, consttype_node, 0, TRUE },
  { "NODE_DYNAMICMODE", NODE_DYNAMICMODE, consttype_node, 0, TRUE },
  { "NODE_RESTARTMODE", NODE_RESTARTMODE, consttype_node, 0, TRUE },
  { "NODE_BREADTHFIRSTMODE", NODE_BREADTHFIRSTMODE, consttype_node, 0, TRUE },
  { "NODE_AUTOORDER", NODE_AUTOORDER, consttype_node, 0, TRUE },
  { "NODE_RCOSTFIXING", NODE_RCOSTFIXING, consttype_node, 0, TRUE },
  { "NODE_STRONGINIT", NODE_STRONGINIT, consttype_node, 0, TRUE },

  { "PRESOLVE_NONE", PRESOLVE_NONE, consttype_presolve, ~0, TRUE },
  { "PRESOLVE_ROWS", PRESOLVE_ROWS, consttype_presolve, 0, TRUE },
  { "PRESOLVE_COLS", PRESOLVE_COLS, consttype_presolve, 0, TRUE },
  { "PRESOLVE_LINDEP", PRESOLVE_LINDEP, consttype_presolve, 0, TRUE },
  { "PRESOLVE_SOS", PRESOLVE_SOS, consttype_presolve, 0, TRUE },
  { "PRESOLVE_REDUCEMIP", PRESOLVE_REDUCEMIP, consttype_presolve, 0, TRUE },
  { "PRESOLVE_KNAPSACK", PRESOLVE_KNAPSACK, consttype_presolve, 0, TRUE },
  { "PRESOLVE_ELIMEQ2", PRESOLVE_ELIMEQ2, consttype_presolve, 0, TRUE },
  { "PRESOLVE_IMPLIEDFREE", PRESOLVE_IMPLIEDFREE, consttype_presolve, 0, TRUE },
  { "PRESOLVE_REDUCEGCD", PRESOLVE_REDUCEGCD, consttype_presolve, 0, TRUE },
  { "PRESOLVE_PROBEFIX", PRESOLVE_PROBEFIX, consttype_presolve, 0, TRUE },
  { "PRESOLVE_PROBEREDUCE", PRESOLVE_PROBEREDUCE, consttype_presolve, 0, TRUE },
  { "PRESOLVE_ROWDOMINATE", PRESOLVE_ROWDOMINATE, consttype_presolve, 0, TRUE },
  { "PRESOLVE_COLDOMINATE", PRESOLVE_COLDOMINATE, consttype_presolve, 0, TRUE },
  { "PRESOLVE_MERGEROWS", PRESOLVE_MERGEROWS, consttype_presolve, 0, TRUE },
  { "PRESOLVE_IMPLIEDSLK", PRESOLVE_IMPLIEDSLK, consttype_presolve, 0, TRUE },
  { "PRESOLVE_COLFIXDUAL", PRESOLVE_COLFIXDUAL, consttype_presolve, 0, TRUE },
  { "PRESOLVE_BOUNDS", PRESOLVE_BOUNDS, consttype_presolve, 0, TRUE },
  { "PRESOLVE_DUALS", PRESOLVE_DUALS, consttype_presolve, 0, TRUE },
  { "PRESOLVE_SENSDUALS", PRESOLVE_SENSDUALS, consttype_presolve, 0, TRUE },

  { "PRICER_FIRSTINDEX", PRICER_FIRSTINDEX, consttype_pricer, 0x03, TRUE },
  { "PRICER_DANTZIG", PRICER_DANTZIG, consttype_pricer, 0x03, TRUE },
  { "PRICER_DEVEX", PRICER_DEVEX, consttype_pricer, 0x03, TRUE },
  { "PRICER_STEEPESTEDGE", PRICER_STEEPESTEDGE, consttype_pricer, 0x03, TRUE },

  { "PRICE_PRIMALFALLBACK", PRICE_PRIMALFALLBACK, consttype_price, 0, TRUE },
  { "PRICE_MULTIPLE", PRICE_MULTIPLE, consttype_price, 0, TRUE },
  { "PRICE_PARTIAL", PRICE_PARTIAL, consttype_price, 0, TRUE },
  { "PRICE_ADAPTIVE", PRICE_ADAPTIVE, consttype_price, 0, TRUE },
  { "PRICE_RANDOMIZE", PRICE_RANDOMIZE, consttype_price, 0, TRUE },
  { "PRICE_AUTOPARTIAL", PRICE_AUTOPARTIAL, consttype_price, 0, TRUE },
  { "PRICE_LOOPLEFT", PRICE_LOOPLEFT, consttype_price, 0, TRUE },
  { "PRICE_LOOPALTERNATE", PRICE_LOOPALTERNATE, consttype_price, 0, TRUE },
  { "PRICE_HARRISTWOPASS", PRICE_HARRISTWOPASS, consttype_price, 0, TRUE },
  { "PRICE_TRUENORMINIT", PRICE_TRUENORMINIT, consttype_price, 0, TRUE },

  { "SCALE_NONE", SCALE_NONE, consttype_scale, 0x07, TRUE },
  { "SCALE_EXTREME", SCALE_EXTREME, consttype_scale, 0x07, TRUE },
  { "SCALE_RANGE", SCALE_RANGE, consttype_scale, 0x07, TRUE },
  { "SCALE_MEAN", SCALE_MEAN, consttype_scale, 0x07, TRUE },
  { "SCALE_GEOMETRIC", SCALE_GEOMETRIC, consttype_scale, 0x07, TRUE },
  { "SCALE_CURTISREID", SCALE_CURTISREID, consttype_scale, 0x07, TRUE },
  { "SCALE_QUADRATIC", SCALE_QUADRATIC, consttype_scale, 0, TRUE },
  { "SCALE_LOGARITHMIC", SCALE_LOGARITHMIC, consttype_scale, 0, TRUE },
  { "SCALE_USERWEIGHT", SCALE_USERWEIGHT, consttype_scale, 0, TRUE },
  { "SCALE_POWER2", SCALE_POWER2, consttype_scale, 0, TRUE },
  { "SCALE_EQUILIBRATE", SCALE_EQUILIBRATE, consttype_scale, 0, TRUE },
  { "SCALE_INTEGERS", SCALE_INTEGERS, consttype_scale, 0, TRUE },
  { "SCALE_DYNUPDATE", SCALE_DYNUPDATE, consttype_scale, 0, TRUE },
  { "SCALE_ROWSONLY", SCALE_ROWSONLY, consttype_scale, 0, TRUE },
  { "SCALE_COLSONLY", SCALE_COLSONLY, consttype_scale, 0, TRUE },

  { "SIMPLEX_PRIMAL_PRIMAL", SIMPLEX_PRIMAL_PRIMAL, consttype_simplex, ~0, TRUE },
  { "SIMPLEX_DUAL_PRIMAL", SIMPLEX_DUAL_PRIMAL, consttype_simplex, ~0, TRUE },
  { "SIMPLEX_PRIMAL_DUAL", SIMPLEX_PRIMAL_DUAL, consttype_simplex, ~0, TRUE },
  { "SIMPLEX_DUAL_DUAL", SIMPLEX_DUAL_DUAL, consttype_simplex, ~0, TRUE },

  { "EPS_TIGHT", EPS_TIGHT, consttype_eps, ~0, TRUE },
  { "EPS_MEDIUM", EPS_MEDIUM, consttype_eps, ~0, TRUE },
  { "EPS_LOOSE", EPS_LOOSE, consttype_eps, ~0, TRUE },
  { "EPS_BAGGY", EPS_BAGGY, consttype_eps, ~0, TRUE },
};

static void createconstant(structlpsolve *lpsolve, int value, unsigned int consttype, char *buf)
{
        int i;

        *buf = 0;
        for (i = 0; i < (int) (sizeof(constants)/sizeof(*constants)); i++)
                if ((constants[i].usereturn) && ((constants[i].consttype & consttype) != 0)) {
                        if (((constants[i].mask == 0) && ((constants[i].value & value) == constants[i].value)) ||
                            ((constants[i].mask != 0) && (constants[i].value == (value & constants[i].mask)))) {
                                if (*buf)
                                        strcat(buf, "|");
                                strcat(buf, constants[i].svalue);
                        }
                }
}

static void returnconstant(structlpsolve *lpsolve, int value, unsigned int consttype)
{
        if (return_constants) {
                char buf[bufsz], *ptr = buf;

                createconstant(lpsolve, value, consttype, buf);
                CreateString(&lpsolve->lpsolvecaller, &ptr, 1, 0);
        }
        else {
                Long    *ipr;

                ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        	*ipr = value;
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
        }
}

static int constantfromstr(structlpsolve *lpsolve, char *buf, unsigned int consttype)
{
        hashelem *hp;
        int value = 0;
        char *ptr, *ptr1, *ptr0 = buf;

        while (*ptr0) {
                for (ptr = ptr0; (*ptr) && (*ptr != '|'); ptr++)
                        *ptr = toupper(*ptr);
                if (*ptr == '|') {
                        *ptr = 0;
                        ptr1 = ptr + 1;
                }
                else
                        ptr1 = ptr;
                while ((*ptr0) && (isspace(*ptr0)))
                        ptr0++;
                ptr--;
                while ((ptr >= ptr0) && (isspace(*ptr)))
                        *(ptr--) = 0;
                hp = findhash(ptr0, constanthash);
                if (hp == NULL) {
                        strcpy(buf, ptr0);
                        strncat(buf,": Unknown.", bufsz);
                        ErrMsgTxt(&lpsolve->lpsolvecaller, buf);
                }
                if ((constants[hp->index].consttype & consttype) == 0) {
                        strcpy(buf, ptr0);
                        strncat(buf,": Not allowed here.", bufsz);
                        ErrMsgTxt(&lpsolve->lpsolvecaller, buf);
                }

                if ((constants[hp->index].mask) && ((value & constants[hp->index].mask) != 0)) {
                        strcpy(buf, ptr0);
                        strncat(buf," cannot be combined with ", bufsz);
                        createconstant(lpsolve, value & constants[hp->index].mask, consttype, buf + strlen(buf));
                        ErrMsgTxt(&lpsolve->lpsolvecaller, buf);
                }

                value |= constants[hp->index].value;

                ptr0 = ptr1;
        }

        return(value);
}

static int constant(structlpsolve *lpsolve, int element, unsigned int consttype)
{
        int value;
        char buf[bufsz];

        if (GetString(&lpsolve->lpsolvecaller, NULL, element, buf, bufsz, FALSE))
                value = constantfromstr(lpsolve, buf, consttype);
        else
                value = (int) GetRealScalar(&lpsolve->lpsolvecaller, element);
        return(value);
}

typedef void (impl_routine)(structlpsolve *lpsolve);

static void impl_set_obj_fn(structlpsolve *lpsolve);

static void addallocmem(structlpsolve *lpsolve, void *ptr)
{
        struct structallocatedmemory *allocatedmemory;

        allocatedmemory = (struct structallocatedmemory *) matCalloc(1, sizeof(*allocatedmemory));
        allocatedmemory->ptr = ptr;
        allocatedmemory->next = lpsolve->allocatedmemory;
        lpsolve->allocatedmemory = allocatedmemory;
}

static void freeallocmem(structlpsolve *lpsolve, void *ptr)
{
        struct structallocatedmemory *allocatedmemory, *allocatedmemory0 = NULL;

        for (allocatedmemory = lpsolve->allocatedmemory; allocatedmemory != NULL; allocatedmemory = allocatedmemory->next) {
                if (allocatedmemory->ptr == ptr) {
                        if (allocatedmemory0 == NULL)
                                lpsolve->allocatedmemory = allocatedmemory->next;
                        else
                                allocatedmemory0->next = allocatedmemory->next;
                        matFree(allocatedmemory);
                        break;
                }
                allocatedmemory0 = allocatedmemory;
        }
}

static void *mallocmem(structlpsolve *lpsolve, size_t size)
{
        void *ptr;

        ptr = matCalloc(size, 1);

        addallocmem(lpsolve, ptr);

        return(ptr);
}

static void *callocmem(structlpsolve *lpsolve, size_t nitems, size_t size)
{
        void *ptr;

        ptr = matCalloc(nitems, size);

        addallocmem(lpsolve, ptr);

        return(ptr);
}

static void *reallocmem(structlpsolve *lpsolve, void *ptr, size_t size)
{
        if (ptr != NULL)
                freeallocmem(lpsolve, ptr);
        ptr = realloc(ptr, size);
        addallocmem(lpsolve, ptr);

        return(ptr);
}

static void freemem(structlpsolve *lpsolve, void *ptr)
{
        if (ptr != NULL) {
                freeallocmem(lpsolve, ptr);
                matFree(ptr);
        }
}

void exitnow(structlpsolvecaller *lpsolvecaller)
{
        longjmp(lpsolvecaller->exit_mark, -1);
}

int EndOfPgr(int x)
{
        return(0);
}

static void Check_nrhs(structlpsolve *lpsolve, int nrhs0)
{
	if (lpsolve->lpsolvecaller.nrhs - 1 != nrhs0) {
                char buf[bufsz];

                sprintf(buf, "%s requires %d argument%s.", lpsolve->cmd, nrhs0, (nrhs0 == 1) ? "" : "s");
		ErrMsgTxt(&lpsolve->lpsolvecaller, buf);
	}
}


/* callback function for lp_solve. Messages are reported via this routine and printed in the application */

static void __WINAPI mylog(lprec *lp, void *userhandle, char *buf)
{
  	Printf("%s", buf);
    mexEvalString("drawnow;"); //flush draw buffer
}


static int __WINAPI myabort(lprec *lp, void *userhandle)
{
        CheckInterrupted((structlpsolvecaller *) userhandle, &interrupted);
        return(interrupted);
}


/* put lp on list */

static int create_handle(structlpsolve *lpsolve, lprec *lp0, char *err)
{
        int i;

        if (lp0 == NULL)
        	ErrMsgTxt(&lpsolve->lpsolvecaller, err);
	for (i = 0; (i <= lp_last) && (lp[i] != NULL); i++);
	if (i > lp_last) {
	  	i = ++lp_last;
                if ((i % LPSTEP) == 0) {
#if defined FORTIFY
                        Fortify_Disable(-1);
#endif
                        if (i == 0)
                                lp = (lprec **) malloc(LPSTEP * sizeof(*lp));
                        else
                                lp = (lprec **) realloc(lp, (i + LPSTEP) * sizeof(*lp));
#if defined FORTIFY
                        Fortify_Disable(-2);
#endif
                        memset(lp + i, 0, LPSTEP * sizeof(*lp));
                }
        }
        lp[i] = lp0;

        putlogfunc(lp0, mylog, &lpsolve->lpsolvecaller);
        set_outputfile(lp0, "");
        putabortfunc(lp0, myabort, &lpsolve->lpsolvecaller);

        return(i);
}

/* lp validation test */

static int handle_valid(int handle)
{
        if (handle < 0 || handle > lp_last || lp[handle] == NULL)
	   return(0);
        else
           return(1);
}

/* free lp from list */

static void delete_handle(int handle)
{
        if (handle_valid(handle)) {
                char *name;
                lprec *lp0 = lp[handle];

                name = get_lp_name(lp0);
                if ((handlehash != NULL) && (name != NULL) && (*name) && (strcmp(name, "Unnamed"))) {
#if defined FORTIFY
                Fortify_Disable(-1);
#endif
                        drophash(name, NULL, handlehash);
#if defined FORTIFY
                Fortify_Disable(-2);
#endif
                }
                delete_lp(lp0);
                lp[handle] = NULL;
        }
}


static void set_handlename(lprec *lp0, char *name, int h)
{
        if (*name) {
#if defined FORTIFY
                Fortify_Disable(-1);
#endif
        	if (handlehash == NULL)
        		handlehash = create_hash_table(100, 0);
                else {
                        char *oldname;

                        oldname = get_lp_name(lp0);
                        if ((oldname != NULL) && (*oldname) && (strcmp(oldname, "Unnamed")))
                                drophash(oldname, NULL, handlehash);
                }
                if (findhash(name, handlehash) == NULL)
                	puthash(name, h, NULL, handlehash);
#if defined FORTIFY
                Fortify_Disable(-2);
#endif
        }
}


/* An exit function to clean up the lp structure.
   called on exit */

void ExitFcn(void)
{
	int	i;

        if (initialized) {
        	for (i = 0; i <= lp_last; i++)
                        delete_handle(i);
#if defined FORTIFY
                Fortify_Disable(-1);
#endif
                free_hash_table(constanthash);
                free_hash_table(cmdhash);
                if (handlehash != NULL)
                	free_hash_table(handlehash);
#if defined FORTIFY
                Fortify_Disable(-2);
#endif
                exit_lpsolve_lib();
#               if defined DEBUG
		        Printf("Terminated%s", NEWLINE);
#               endif
        }
}


#if defined DEBUG || defined DEMO
/* xxlpsolve('demo') */
/* This demo is not a necessary part of the program */
#endif

static void impl_demo(structlpsolve *lpsolve)
{
        int i, h;
        lprec *lp0;

	/* Printf("%.15f", GetRealScalar(&lpsolve->lpsolvecaller, 1)); */

        Check_nrhs(lpsolve, 0);

        h = create_handle(lpsolve, make_lp(0, 4), "make_lp failed");
        lp0 = lp[h];

        Printf("min: 2 C1 +3 C2 -2 C3 +3 C4;%s", NEWLINE);
	str_set_obj_fn(lp0, "2 3 -2 3");

        Printf("3 C1 +2 C2 +2 C3 +1 C4 <= 4.0;%s", NEWLINE);
        str_add_constraint(lp0, "3 2 2 1", LE, 4.0);

        Printf("0 C1 +4 C2 +3 C3 +1 C4 >= 3.0;%s", NEWLINE);
	str_add_constraint(lp0, "0 4 3 1", GE, 3.0);

        set_verbose(lp0, 0);
        Printf("solve: %d%s", solve(lp0), NEWLINE);
        Printf("obj: %f%s", get_objective(lp0), NEWLINE);
        for (i = 1; i <= get_Ncolumns(lp0); i++)
        	Printf("C%d: %f%s", i, get_var_primalresult(lp0, get_Nrows(lp0) + i), NEWLINE);
        write_lp(lp0, "a.lp");

        delete_handle(h);
}


/* since always a matrix vector is provided to both add_column and add_columnex, always call the more
   performant sparse version of the two routines */
/* return = xxlpsolve('add_column', lp, [column]) */
/* return = xxlpsolve('add_columnex', lp, [column]) */

static void impl_add_column(structlpsolve *lpsolve)
{
        int m, count;
        int	*index;
	REAL	*vec;
        int     result;
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        m = get_Nrows(lpsolve->lp);
        vec = (REAL *) callocmem(lpsolve, 1 + m, sizeof(*vec));
        index = (int *) callocmem(lpsolve, 1 + m, sizeof(*index));
	count = GetRealSparseVector(&lpsolve->lpsolvecaller, 2, vec, index, 0, 1 + m, 0);
	result = add_columnex(lpsolve->lp, count, vec, index);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
        freemem(lpsolve, index);
        freemem(lpsolve, vec);
}


/* since always a matrix vector is provided to both add_constraint and add_constraintex, always call the more
   performant sparse version of the two routines */
/* return = xxlpsolve('add_constraint', lp, [row], constr_type, rh) */
/* return = xxlpsolve('add_constraintex', lp, [row], constr_type, rh) */

static void impl_add_constraint(structlpsolve *lpsolve)
{
        int type, n, count;
        int	*index;
	REAL	*vec, value;
        int     result;
        Long    *ipr;

        Check_nrhs(lpsolve, 4);

        type = constant(lpsolve, 3, consttype_constrainttype);
        value = GetRealScalar(&lpsolve->lpsolvecaller, 4);
        n = get_Ncolumns(lpsolve->lp);
        vec = (REAL *) callocmem(lpsolve, n, sizeof(*vec));
        index = (int *) callocmem(lpsolve, n, sizeof(*index));
	count = GetRealSparseVector(&lpsolve->lpsolvecaller, 2, vec, index, 1, n, 0);
	result = add_constraintex(lpsolve->lp, count, vec, index, type, value);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
        freemem(lpsolve, index);
        freemem(lpsolve, vec);
}


/* return = xxlpsolve('add_SOS', lp, name, sostype, priority, [sosvars], [weights]) */

static void impl_add_SOS(structlpsolve *lpsolve)
{
        int n, *sosvars;
        REAL *weights;
        int count1, count2;
        int     result;
        char buf[bufsz];
        Long    *ipr;

        Check_nrhs(lpsolve, 6);
        GetString(&lpsolve->lpsolvecaller, NULL, 2, buf, bufsz, TRUE);
        n = get_Ncolumns(lpsolve->lp);
        sosvars = (int *) callocmem(lpsolve, n, sizeof(*sosvars));
        weights = (REAL *) callocmem(lpsolve, n, sizeof(*weights));
	count1 = GetIntVector(&lpsolve->lpsolvecaller, 5, sosvars, 0, n, FALSE);
        count2 = GetRealVector(&lpsolve->lpsolvecaller, 6, weights, 0, n, FALSE);
        if (count1 != count2) {
          freemem(lpsolve, weights);
          freemem(lpsolve, sosvars);
          ErrMsgTxt(&lpsolve->lpsolvecaller, "add_SOS: sosvars and weights vector must have same size.");
        }
	result = add_SOS(lpsolve->lp, buf, (int) GetRealScalar(&lpsolve->lpsolvecaller, 3), (int) GetRealScalar(&lpsolve->lpsolvecaller, 4), count1, sosvars,weights);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
        freemem(lpsolve, weights);
        freemem(lpsolve, sosvars);
}


/* return = xxlpsolve('column_in_lp', lp, [column]) */

static void impl_column_in_lp(structlpsolve *lpsolve)
{
        int n;
	REAL	*vec;
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        n = get_Nrows(lpsolve->lp);
        vec = (REAL *) callocmem(lpsolve, 1 + n, sizeof(*vec));
        GetRealVector(&lpsolve->lpsolvecaller, 2, vec, 0, 1 + n, TRUE);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = column_in_lp(lpsolve->lp, vec);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
        freemem(lpsolve, vec);
}


/* xxlpsolve('copy_lp', lp) */

static void impl_copy_lp(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = create_handle(lpsolve, copy_lp(lpsolve->lp), "copy_lp failed");
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* xxlpsolve('default_basis', lp) */

static void impl_default_basis(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 1);
        default_basis(lpsolve->lp);
}


/* return = xxlpsolve('del_column', lp, column) */

static void impl_del_column(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = del_column(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2));
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('del_constraint', lp, del_row) */

static void impl_del_constraint(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = del_constraint(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2));
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* xxlpsolve('delete_lp', lp) */
/* xxlpsolve('free_lp', lp) */

static void impl_delete_lp(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 1);
        delete_handle(lpsolve->h);
}


/* xxlpsolve('dualize_lp', lp) */

static void impl_dualize_lp(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = dualize_lp(lpsolve->lp);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_anti_degen', lp) */
static void impl_get_anti_degen(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 1);
        returnconstant(lpsolve, get_anti_degen(lpsolve->lp), consttype_antidegen);
}


/* [bascolumn] = xxlpsolve('get_basis', lp {, nonbasic}) */

static void impl_get_basis(structlpsolve *lpsolve)
{
        int n, i, *bascolumn, *bascolumn0;
        MYBOOL nonbasic;
        Long    *ipr, *ipr0;

        if (lpsolve->lpsolvecaller.nrhs == 1+1)
                n = 1;
        else
                n = 2;
        Check_nrhs(lpsolve, n);
        if (n == 1)
                nonbasic = 0;
        else
        	nonbasic = (MYBOOL) GetRealScalar(&lpsolve->lpsolvecaller, 2);
        n = get_Nrows(lpsolve->lp) + ((nonbasic) ? get_Ncolumns(lpsolve->lp) : 0);
        bascolumn0 = bascolumn = (int *) callocmem(lpsolve, 1 + n, sizeof(*bascolumn));
        if(get_basis(lpsolve->lp, bascolumn, nonbasic)) {
                ipr0 = ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, n, 1, 0);
                for (i = 0; i < n; i++)
                  *(ipr++) = *(++bascolumn);
        }
        else {
                n = 0;
                ipr0 = ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, n, 1, 0);
        }
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr0, n, 1, 0, TRUE);
        freemem(lpsolve, bascolumn0);
}


/* return = xxlpsolve('get_basiscrash', lp) */

static void impl_get_basiscrash(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 1);
        returnconstant(lpsolve, get_basiscrash(lpsolve->lp), consttype_crash);
}


/* return = xxlpsolve('get_bb_depthlimit', lp) */

static void impl_get_bb_depthlimit(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = get_bb_depthlimit(lpsolve->lp);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_bb_floorfirst', lp) */

static void impl_get_bb_floorfirst(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 1);
        returnconstant(lpsolve, get_bb_floorfirst(lpsolve->lp), consttype_branch);
}


/* return = xxlpsolve('get_bb_rule', lp) */

static void impl_get_bb_rule(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 1);
        returnconstant(lpsolve, get_bb_rule(lpsolve->lp), consttype_node);
}


/* return = xxlpsolve('get_bounds_tighter', lp) */

static void impl_get_bounds_tighter(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = get_bounds_tighter(lpsolve->lp);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_break_at_value', lp) */

static void impl_get_break_at_value(structlpsolve *lpsolve)
{
        Double  *dpr;

        Check_nrhs(lpsolve, 1);
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *dpr = get_break_at_value(lpsolve->lp);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, 1, 0, TRUE);
}


/* name = xxlpsolve('get_col_name', lp, column) */
/* [names] = xxlpsolve('get_col_name', lp) */

static void impl_get_col_name(structlpsolve *lpsolve)
{
        char *name;

        if (lpsolve->lpsolvecaller.nrhs == 1+1) {
                int n, i;
                char **names;

                n = get_Ncolumns(lpsolve->lp);

                names = (char **) callocmem(lpsolve, n, sizeof(*names));
                for (i = 0; i < n; i++) {
                        name = get_col_name(lpsolve->lp, i + 1);
                        if (name == NULL)
                                name = "";
                        names[i] = (char *) mallocmem(lpsolve, strlen(name) + 1);
                        strcpy(names[i], name);
                }
                CreateString(&lpsolve->lpsolvecaller, names, n, 0);
                for (i = 0; i < n; i++)
                        freemem(lpsolve, names[i]);
                freemem(lpsolve, names);
        }
        else {
        	Check_nrhs(lpsolve, 2);
                name = get_col_name(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2));
                if (name == NULL)
                        name = "";
                CreateString(&lpsolve->lpsolvecaller, &name, 1, 0);
        }
}


/* [column, return] = xxlpsolve('get_column', lp, col_nr) */
/* [column, return] = xxlpsolve('get_columnex', lp, col_nr) */

static void impl_get_column(structlpsolve *lpsolve)
{
        int col;
        int     result;
        Double  *dpr;
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
	col = (int) GetRealScalar(&lpsolve->lpsolvecaller, 2);
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1 + get_Nrows(lpsolve->lp), 1, 0);
	result = get_column(lpsolve->lp, col, dpr);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1 + get_Nrows(lpsolve->lp), 1, 0, TRUE);
        if (lpsolve->lpsolvecaller.nlhs > 1) {
                ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 1);
                *ipr = result;
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 1, TRUE);
        }
}


/* return = xxlpsolve('get_constr_type', lp, row) */
/* [constr_type] = xxlpsolve('get_constr_type', lp) */

static void impl_get_constr_type(structlpsolve *lpsolve)
{
        Long    *ipr, *ipr0;

	if (lpsolve->lpsolvecaller.nrhs == 1+1) {
                int i, m;

                Check_nrhs(lpsolve, 1);
                m = get_Nrows(lpsolve->lp);
                if (return_constants) {
                	char **types = NULL, *ptr, buf[bufsz];

                        types = (char **) callocmem(lpsolve, m, sizeof(*types));
                        for (i = 1; i <= m; i++) {
                                createconstant(lpsolve, get_constr_type(lpsolve->lp, i), consttype_constrainttype, buf);
                                ptr = types[i - 1] = (char *) callocmem(lpsolve, strlen(buf) + 1, sizeof(**types));
                                strcpy(ptr, buf);
                        }
                        CreateString(&lpsolve->lpsolvecaller, types, m, 0);
                        for (i = 0; i < m; i++)
                                freemem(lpsolve, types[i]);
                        freemem(lpsolve, types);
                }
                else {
                        ipr0 = ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, m, 1, 0);
                        for (i = 1; i <= m; i++)
                                *(ipr++) = get_constr_type(lpsolve->lp, i);
                        SetLongMatrix(&lpsolve->lpsolvecaller, ipr0, m, 1, 0, TRUE);
                }
        }
        else {
        	Check_nrhs(lpsolve, 2);
                returnconstant(lpsolve, get_constr_type(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2)), consttype_constrainttype);
        }
}

/* return = xxlpsolve('get_constr_value', lp, row, [primsolution]) */

static void impl_get_constr_value(structlpsolve *lpsolve)
{
        int n, count;
        int	*nzindex;
	REAL	*primsolution;
        Double  *dpr;

        if (lpsolve->lpsolvecaller.nrhs == 1+2) {
                Check_nrhs(lpsolve, 2);
                primsolution = NULL;
                nzindex = NULL;
                count = 0;
        }
        else {
                Check_nrhs(lpsolve, 3);
                n = get_Ncolumns(lpsolve->lp);
                if (n == 0)
                        n = 1;
                primsolution = (REAL *) callocmem(lpsolve, n, sizeof(*primsolution));
                nzindex = (int *) callocmem(lpsolve, n, sizeof(*nzindex));
        	count = GetRealSparseVector(&lpsolve->lpsolvecaller, 3, primsolution, nzindex, 1, n, 0);
        }
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *dpr = get_constr_value(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2), count, primsolution, nzindex);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, 1, 0, TRUE);
        if (nzindex != NULL)
                freemem(lpsolve, nzindex);
        if (primsolution != NULL)
                freemem(lpsolve, primsolution);
}


/* [constr, return] = xxlpsolve('get_constraints', lp) */

static void impl_get_constraints(structlpsolve *lpsolve)
{
        int     result;
        Double  *dpr;
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, get_Nrows(lpsolve->lp), 1, 0);
	result = get_constraints(lpsolve->lp, dpr);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, get_Nrows(lpsolve->lp), 1, 0, TRUE);
        if (lpsolve->lpsolvecaller.nlhs > 1) {
                ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 1);
                *ipr = result;
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 1, TRUE);
        }
}


/* [duals, return] = xxlpsolve('get_dual_solution', lp) */

static void impl_get_dual_solution(structlpsolve *lpsolve)
{
        int i;
	REAL	*vec = NULL;
        int     result;
        Double  *dpr;
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
	result = get_ptr_dual_solution(lpsolve->lp, &vec);
        if ((!result) || (vec == NULL))
                ErrMsgTxt(&lpsolve->lpsolvecaller, "get_dual_solution: sensitivity unknown.");
        i = get_Nrows(lpsolve->lp) + get_Ncolumns(lpsolve->lp);
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, i, 1, 0);
        memcpy(dpr, vec + 1, i * sizeof(*vec));
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, i, 1, 0, TRUE);
        if (lpsolve->lpsolvecaller.nlhs > 1) {
                ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 1);
                *ipr = result;
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 1, TRUE);
        }
}


/* return = xxlpsolve('get_epsb', lp) */

static void impl_get_epsb(structlpsolve *lpsolve)
{
        Double  *dpr;

        Check_nrhs(lpsolve, 1);
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *dpr = get_epsb(lpsolve->lp);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_epsd', lp) */

static void impl_get_epsd(structlpsolve *lpsolve)
{
        Double  *dpr;

        Check_nrhs(lpsolve, 1);
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *dpr = get_epsd(lpsolve->lp);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_epsel', lp) */

static void impl_get_epsel(structlpsolve *lpsolve)
{
        Double  *dpr;

        Check_nrhs(lpsolve, 1);
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *dpr = get_epsel(lpsolve->lp);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_epsint', lp) */

static void impl_get_epsint(structlpsolve *lpsolve)
{
        Double  *dpr;

        Check_nrhs(lpsolve, 1);
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *dpr = get_epsint(lpsolve->lp);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_epsperturb', lp) */

static void impl_get_epsperturb(structlpsolve *lpsolve)
{
        Double  *dpr;

        Check_nrhs(lpsolve, 1);
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *dpr = get_epsperturb(lpsolve->lp);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_epspivot', lp) */

static void impl_get_epspivot(structlpsolve *lpsolve)
{
        Double  *dpr;

        Check_nrhs(lpsolve, 1);
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *dpr = get_epspivot(lpsolve->lp);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_improve', lp) */

static void impl_get_improve(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 1);
        returnconstant(lpsolve, get_improve(lpsolve->lp), consttype_improve);
}


/* return = xxlpsolve('get_infinite', lp) */

static void impl_get_infinite(structlpsolve *lpsolve)
{
        Double  *dpr;

        Check_nrhs(lpsolve, 1);
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *dpr = get_infinite(lpsolve->lp);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_lowbo', lp, column) */
/* [return] = xxlpsolve('get_lowbo', lp) */

static void impl_get_lowbo(structlpsolve *lpsolve)
{
        Double  *dpr, *dpr0;

	if (lpsolve->lpsolvecaller.nrhs == 1+1) {
                int n, i;

                Check_nrhs(lpsolve, 1);
                n = get_Ncolumns(lpsolve->lp);
                dpr0 = dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, n, 1, 0);
                for (i = 1; i <= n; i++)
                  *(dpr++) = get_lowbo(lpsolve->lp, i);
                SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr0, n, 1, 0, TRUE);
        }
        else {
        	Check_nrhs(lpsolve, 2);
                dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
		*dpr = get_lowbo(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2));
                SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, 1, 0, TRUE);
        }
}


/* return = xxlpsolve('get_lp_index', lp, orig_index) */

static void impl_get_lp_index(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = get_lp_index(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2));
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* name = xxlpsolve('get_lp_name', lp) */

static void impl_get_lp_name(structlpsolve *lpsolve)
{
        char *name;

        Check_nrhs(lpsolve, 1);
        name = get_lp_name(lpsolve->lp);
        CreateString(&lpsolve->lpsolvecaller, &name, 1, 0);
}


/* value = xxlpsolve('get_mat', lp, row, col) */
/* [matrix, return] = xxlpsolve('get_mat', lp[, sparse]) */

static void impl_get_mat(structlpsolve *lpsolve)
{
        Double  *dpr, *dpr0;
        Long    *ipr;

	if ((lpsolve->lpsolvecaller.nrhs == 1+1) || (lpsolve->lpsolvecaller.nrhs == 1+2)) {
                int m, n, i, sparse = FALSE;
        	REAL	*vec;
                int     result;

                Check_nrhs(lpsolve, lpsolve->lpsolvecaller.nrhs - 1);
                m = get_Nrows(lpsolve->lp);
                n = get_Ncolumns(lpsolve->lp);
                vec = (REAL *) callocmem(lpsolve, 1 + m, sizeof(*vec));
                result = TRUE;
                if (lpsolve->lpsolvecaller.nrhs == 1+2)
                        sparse = (int) GetRealScalar(&lpsolve->lpsolvecaller, 2);
                if (!sparse) { /* return the matrix in dense format */
                        dpr0 = dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, m, n, 0);
                        for (i = 1; (i <= n) && (result); i++) {
                		result = get_column(lpsolve->lp, i, vec);
                        	memcpy(dpr, vec + 1, m * sizeof(*vec));
                                dpr += m;
                        }
                }
                else {         /* return the matrix in sparse format */
                        int nz = 0;

                        dpr0 = dpr = CreateDoubleSparseMatrix(&lpsolve->lpsolvecaller, m, n, 0);
                        for (i = 1; (i <= n) && (result); i++) {
                		result = get_column(lpsolve->lp, i, vec);
                                SetColumnDoubleSparseMatrix(&lpsolve->lpsolvecaller, 0, m, n, dpr, i, vec + 1, NULL, m, &nz);
                        }
                }
                SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr0, m, n, 0, TRUE);
                freemem(lpsolve, vec);
                if (lpsolve->lpsolvecaller.nlhs > 1) {
                        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 1);
                        *ipr = result;
                        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 1, TRUE);
                }
        }
        else {
	        Check_nrhs(lpsolve, 3);
                dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        	*dpr = get_mat(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2), (int) GetRealScalar(&lpsolve->lpsolvecaller, 3));
                SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, 1, 0, TRUE);
        }
}


/* return = xxlpsolve('get_max_level', lp) */

static void impl_get_max_level(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = get_max_level(lpsolve->lp);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_maxpivot', lp) */

static void impl_get_maxpivot(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = get_maxpivot(lpsolve->lp);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_mip_gap', lp, absolute) */

static void impl_get_mip_gap(structlpsolve *lpsolve)
{
        Double  *dpr;

        Check_nrhs(lpsolve, 2);
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *dpr = get_mip_gap(lpsolve->lp, (MYBOOL) GetRealScalar(&lpsolve->lpsolvecaller, 2));
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_nameindex', lp, name, isrow) */

static void impl_get_nameindex(structlpsolve *lpsolve)
{
        int     result;
        char buf[bufsz];
        Long    *ipr;

        Check_nrhs(lpsolve, 3);
        GetString(&lpsolve->lpsolvecaller, NULL, 2, buf, bufsz, TRUE);
        result = get_nameindex(lpsolve->lp, buf, (MYBOOL) GetRealScalar(&lpsolve->lpsolvecaller, 3));
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_Ncolumns', lp) */

static void impl_get_Ncolumns(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = get_Ncolumns(lpsolve->lp);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_negrange', lp) */

static void impl_get_negrange(structlpsolve *lpsolve)
{
        Double  *dpr;

        Check_nrhs(lpsolve, 1);
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *dpr = get_negrange(lpsolve->lp);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_nonzeros', lp) */

static void impl_get_nonzeros(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = get_nonzeros(lpsolve->lp);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_Norig_columns', lp) */

static void impl_get_Norig_columns(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = get_Norig_columns(lpsolve->lp);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_Norig_rows', lp) */

static void impl_get_Norig_rows(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = get_Norig_rows(lpsolve->lp);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_Nrows', lp) */

static void impl_get_Nrows(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = get_Nrows(lpsolve->lp);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_obj_bound', lp) */

static void impl_get_obj_bound(structlpsolve *lpsolve)
{
        Double  *dpr;

        Check_nrhs(lpsolve, 1);
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *dpr = get_obj_bound(lpsolve->lp);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, 1, 0, TRUE);
}


/* [row_vec, return] = xxlpsolve('get_obj_fn', lp) */
/* [row_vec, return] = xxlpsolve('get_obj_fun', lp) */

static void impl_get_obj_fn(structlpsolve *lpsolve)
{
        int n;
	REAL	*vec;
        int     result;
        Double  *dpr;
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        n = get_Ncolumns(lpsolve->lp);
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, n, 0);
        vec = (REAL *) callocmem(lpsolve, 1 + n, sizeof(*vec));
	result = get_row(lpsolve->lp, 0, vec);
        memcpy(dpr, vec + 1, n * sizeof(*vec));
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, n, 0, TRUE);
        freemem(lpsolve, vec);
        if (lpsolve->lpsolvecaller.nlhs > 1) {
                ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 1);
                *ipr = result;
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 1, TRUE);
        }
}


/* return = xxlpsolve('get_objective', lp) */

static void impl_get_objective(structlpsolve *lpsolve)
{
        Double  *dpr;

        Check_nrhs(lpsolve, 1);
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *dpr = get_objective(lpsolve->lp);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, 1, 0, TRUE);
}


/* name = xxlpsolve('get_objective_name', lp) */

static void impl_get_objective_name(structlpsolve *lpsolve)
{
        char *name;

        Check_nrhs(lpsolve, 1);
        name = get_row_name(lpsolve->lp, 0);
        CreateString(&lpsolve->lpsolvecaller, &name, 1, 0);
}


/* return = xxlpsolve('get_orig_index', lp, lp_index) */

static void impl_get_orig_index(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = get_orig_index(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2));
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* name = xxlpsolve('get_origcol_name', lp, column) */
/* [names] = xxlpsolve('get_origcol_name', lp) */

static void impl_get_origcol_name(structlpsolve *lpsolve)
{
        char *name;

	if (lpsolve->lpsolvecaller.nrhs == 1+1) {
                int n, i;
        	char **names;

                n = get_Ncolumns(lpsolve->lp);

                names = (char **) callocmem(lpsolve, n, sizeof(*names));
                for (i = 0; i < n; i++) {
                        name = get_origcol_name(lpsolve->lp, i + 1);
                        if (name == NULL)
                                name = "";
                        names[i] = (char *) mallocmem(lpsolve, strlen(name) + 1);
                        strcpy(names[i], name);
                }
                CreateString(&lpsolve->lpsolvecaller, names, n, 0);
                for (i = 0; i < n; i++)
                        freemem(lpsolve, names[i]);
                freemem(lpsolve, names);
        }
        else {
	        Check_nrhs(lpsolve, 2);
                name = get_origcol_name(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2));
                if (name == NULL)
                        name = "";
                CreateString(&lpsolve->lpsolvecaller, &name, 1, 0);
        }
}


/* name = xxlpsolve('get_origrow_name', lp, row) */
/* [names] = xxlpsolve('get_origrow_name', lp) */

static void impl_get_origrow_name(structlpsolve *lpsolve)
{
        char *name;

	if (lpsolve->lpsolvecaller.nrhs == 1+1) {
                int m, i;
        	char **names;

                m = get_Nrows(lpsolve->lp);

                names = (char **) callocmem(lpsolve, m, sizeof(*names));
                for (i = 0; i < m; i++) {
                        name = get_origrow_name(lpsolve->lp, i + 1);
                        if (name == NULL)
                                name = "";
                        names[i] = (char *) mallocmem(lpsolve, strlen(name) + 1);
                        strcpy(names[i], name);
                }
                CreateString(&lpsolve->lpsolvecaller, names, m, 0);
                for (i = 0; i < m; i++)
                        free(names[i]);
                freemem(lpsolve, names);
        }
        else {
        	Check_nrhs(lpsolve, 2);
                name = get_origrow_name(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2));
                if (name == NULL)
                        name = "";
                CreateString(&lpsolve->lpsolvecaller, &name, 1, 0);
        }
}


/* return = xxlpsolve('get_pivoting', lp) */

static void impl_get_pivoting(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 1);
        returnconstant(lpsolve, get_pivoting(lpsolve->lp), consttype_pricer|consttype_price);
}


/* return = xxlpsolve('get_presolve', lp) */

static void impl_get_presolve(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 1);
        returnconstant(lpsolve, get_presolve(lpsolve->lp), consttype_presolve);
}


/* return = xxlpsolve('get_presolveloops', lp) */

static void impl_get_presolveloops(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = get_presolveloops(lpsolve->lp);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* [pv, return] = xxlpsolve('get_primal_solution', lp) */

static void impl_get_primal_solution(structlpsolve *lpsolve)
{
        int i;
        int     result;
        Double  *dpr;
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        i = 1 + get_Nrows(lpsolve->lp) + get_Ncolumns(lpsolve->lp);
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, i, 1, 0);
	result = get_primal_solution(lpsolve->lp, dpr);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, i, 1, 0, TRUE);
        if (lpsolve->lpsolvecaller.nlhs > 1) {
                ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 1);
                *ipr = result;
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 1, TRUE);
        }
}


/* return = xxlpsolve('get_print_sol', lp) */

static void impl_get_print_sol(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = get_print_sol(lpsolve->lp);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_rh', lp, row) */
/* [rh] = xxlpsolve('get_rh', lp) */

static void impl_get_rh(structlpsolve *lpsolve)
{
        Double  *dpr, *dpr0;

	if (lpsolve->lpsolvecaller.nrhs == 1+1) {
                int m, i;

                Check_nrhs(lpsolve, 1);
                m = get_Nrows(lpsolve->lp);
                dpr0 = dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1 + m, 1, 0);
                for (i = 0; i <= m; i++)
                  *(dpr++) = get_rh(lpsolve->lp, i);
                SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr0, 1 + m, 1, 0, TRUE);
        }
        else {
                Check_nrhs(lpsolve, 2);
                dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        	*dpr = get_rh(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2));
                SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, 1, 0, TRUE);
        }
}


/* return = xxlpsolve('get_rh_range', lp, row) */
/* [rh_ranges] = xxlpsolve('get_rh_range', lp) */

static void impl_get_rh_range(structlpsolve *lpsolve)
{
        Double  *dpr, *dpr0;

	if (lpsolve->lpsolvecaller.nrhs == 1+1) {
                int m, i;

                Check_nrhs(lpsolve, 1);
                m = get_Nrows(lpsolve->lp);
                dpr0 = dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, m, 1, 0);
                for (i = 1; i <= m; i++)
                  *(dpr++) = get_rh_range(lpsolve->lp, i);
                SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr0, m, 1, 0, TRUE);
        }
        else {
	        Check_nrhs(lpsolve, 2);
                dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
		*dpr = get_rh_range(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2));
                SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, 1, 0, TRUE);
        }
}


/* [row, return] = xxlpsolve('get_row', lp, row_nr) */
/* [row, return] = xxlpsolve('get_rowex', lp, row_nr) */

static void impl_get_row(structlpsolve *lpsolve)
{
        int n;
	REAL	*vec;
        int     result;
        Double  *dpr;
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        n = get_Ncolumns(lpsolve->lp);
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, n, 0);
        vec = (REAL *) callocmem(lpsolve, 1 + n, sizeof(*vec));
	result = get_row(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2), vec);
        memcpy(dpr, vec + 1, n * sizeof(*vec));
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, n, 0, TRUE);
        freemem(lpsolve, vec);
        if (lpsolve->lpsolvecaller.nlhs > 1) {
                ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 1);
                *ipr = result;
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 1, TRUE);
        }
}


/* name = xxlpsolve('get_row_name', lp, row) */
/* [names] = xxlpsolve('get_row_name', lp) */

static void impl_get_row_name(structlpsolve *lpsolve)
{
        char *name;

        if (lpsolve->lpsolvecaller.nrhs == 1+1) {
                int m, i;
        	char **names;

                m = get_Nrows(lpsolve->lp);

                names = (char **) callocmem(lpsolve, m, sizeof(*names));
                for (i = 0; i < m; i++) {
                        name = get_row_name(lpsolve->lp, i + 1);
                        if (name == NULL)
                                name = "";
                        names[i] = (char *) mallocmem(lpsolve, strlen(name) + 1);
                        strcpy(names[i], name);
                }
                CreateString(&lpsolve->lpsolvecaller, names, m, 0);
                for (i = 0; i < m; i++)
                        freemem(lpsolve, names[i]);
                freemem(lpsolve, names);
        }
        else {
                Check_nrhs(lpsolve, 2);
                name = get_row_name(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2));
                if (name == NULL)
                        name = "";
                CreateString(&lpsolve->lpsolvecaller, &name, 1, 0);
        }
}


/* return = xxlpsolve('get_scalelimit', lp) */

static void impl_get_scalelimit(structlpsolve *lpsolve)
{
        Double  *dpr;

        Check_nrhs(lpsolve, 1);
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *dpr = get_scalelimit(lpsolve->lp);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_scaling', lp) */

static void impl_get_scaling(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 1);
        returnconstant(lpsolve, get_scaling(lpsolve->lp), consttype_scale);
}


/* [objfrom, objtill, objfromvalue, objtillvalue, return] = xxlpsolve('get_sensitivity_obj', lp) */
/* [objfrom, objtill, objfromvalue, objtillvalue, return] = xxlpsolve('get_sensitivity_objex', lp) */

static void impl_get_sensitivity_objex(structlpsolve *lpsolve)
{
        int n;
	REAL	*objfrom = NULL, *objtill = NULL, *objfromvalue, *objtillvalue;
        int     result;
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
	result = get_ptr_sensitivity_obj(lpsolve->lp, &objfrom, &objtill);
        if ((!result) || (objfrom == NULL) || (objtill == NULL))
                ErrMsgTxt(&lpsolve->lpsolvecaller, "get_sensitivity_obj: sensitivity unknown.");
        n = get_Ncolumns(lpsolve->lp);
        objfrom = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, n, 0);
        if (lpsolve->lpsolvecaller.nlhs > 1)
        	objtill = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, n, 1);
        else
                objtill = NULL;
        if (lpsolve->lpsolvecaller.nlhs > 2)
        	objfromvalue = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, n, 2);
        else
                objfromvalue = NULL;
        if (lpsolve->lpsolvecaller.nlhs > 3) {
        	objtillvalue = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, n, 3);
                memset(objtillvalue, 0, n * sizeof(*objtillvalue));
        }
        else
                objtillvalue = NULL;
	result = get_sensitivity_objex(lpsolve->lp, objfrom, objtill, objfromvalue, NULL /* objtillvalue */);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, objfrom, 1, n, 0, TRUE);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, objtill, 1, n, 1, TRUE);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, objfromvalue, 1, n, 2, TRUE);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, objtillvalue, 1, n, 3, TRUE);
        if (lpsolve->lpsolvecaller.nlhs > 4) {
                ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 4);
                *ipr = result;
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 4, TRUE);
        }
}


/* [duals, dualsfrom, dualstill, return] = xxlpsolve('get_sensitivity_rhs', lp) */
/* [duals, dualsfrom, dualstill, return] = xxlpsolve('get_sensitivity_rhsex', lp) */

static void impl_get_sensitivity_rhsex(structlpsolve *lpsolve)
{
        int i;
        REAL *duals = NULL, *dualsfrom = NULL, *dualstill = NULL;
        int     result;
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        result = get_ptr_sensitivity_rhs(lpsolve->lp, &duals, &dualsfrom, &dualstill);
        if ((!result) || (duals == NULL) || (dualsfrom == NULL) || (dualstill == NULL))
                ErrMsgTxt(&lpsolve->lpsolvecaller, "get_sensitivity_rhs: sensitivity unknown.");
        i = get_Nrows(lpsolve->lp) + get_Ncolumns(lpsolve->lp);
        duals = CreateDoubleMatrix(&lpsolve->lpsolvecaller, i, 1, 0);
        if (lpsolve->lpsolvecaller.nlhs > 1)
        	dualsfrom = CreateDoubleMatrix(&lpsolve->lpsolvecaller, i, 1, 1);
        else
		dualsfrom = NULL;
        if (lpsolve->lpsolvecaller.nlhs > 2)
        	dualstill = CreateDoubleMatrix(&lpsolve->lpsolvecaller, i, 1, 2);
        else
                dualstill = NULL;
	result = get_sensitivity_rhs(lpsolve->lp, duals, dualsfrom, dualstill);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, duals, i, 1, 0, TRUE);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dualsfrom, i, 1, 1, TRUE);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dualstill, i, 1, 2, TRUE);
        if (lpsolve->lpsolvecaller.nlhs > 3) {
                ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 3);
                *ipr = result;
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 3, TRUE);
        }
}


/* return = xxlpsolve('get_simplextype', lp) */

static void impl_get_simplextype(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 1);
        returnconstant(lpsolve, get_simplextype(lpsolve->lp), consttype_simplex);
}


/* [obj, x, duals, return] = xxlpsolve('get_solution', lp) */

static void impl_get_solution(structlpsolve *lpsolve)
{
        int m, n;
        REAL *duals;
        int     result = 0;
        Double  *dpr;
        Long    *ipr;

        Check_nrhs(lpsolve, 1);

        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *dpr = get_objective(lpsolve->lp);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, 1, 0, TRUE);

        if (lpsolve->lpsolvecaller.nlhs > 1) {
	        n = get_Ncolumns(lpsolve->lp);
        	dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, n, 1, 1);
	        result = get_variables(lpsolve->lp, dpr);
                SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, n, 1, 1, TRUE);
        }

        if (lpsolve->lpsolvecaller.nlhs > 2) {
                m = get_Nrows(lpsolve->lp);
        	dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, m, 1, 2);
		result &= get_ptr_dual_solution(lpsolve->lp, &duals);
                memcpy(dpr, duals + 1, m * sizeof(*dpr));
                SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, m, 1, 2, TRUE);
        }

        if (lpsolve->lpsolvecaller.nlhs > 3) {
                ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 3);
                *ipr = result;
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 3, TRUE);
        }
}


/* return = xxlpsolve('get_solutioncount', lp) */

static void impl_get_solutioncount(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = get_solutioncount(lpsolve->lp);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_solutionlimit', lp) */

static void impl_get_solutionlimit(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = get_solutionlimit(lpsolve->lp);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_status', lp) */

static void impl_get_status(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = get_status(lpsolve->lp);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_statustext', lp, statuscode) */

static void impl_get_statustext(structlpsolve *lpsolve)
{
        char *name;

        Check_nrhs(lpsolve, 2);
        name = get_statustext(lpsolve->lp, constant(lpsolve, 2, consttype_solve));
	CreateString(&lpsolve->lpsolvecaller, &name, 1, 0);
}


/* return = xxlpsolve('get_timeout', lp) */

static void impl_get_timeout(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = get_timeout(lpsolve->lp);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_total_iter', lp) */

static void impl_get_total_iter(structlpsolve *lpsolve)
{
        Double  *dpr;

        Check_nrhs(lpsolve, 1);
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *dpr = (double) get_total_iter(lpsolve->lp);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_total_nodes', lp) */

static void impl_get_total_nodes(structlpsolve *lpsolve)
{
        Double  *dpr;

        Check_nrhs(lpsolve, 1);
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *dpr = (double) get_total_nodes(lpsolve->lp);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_upbo', lp, column) */
/* [upbo] = xxlpsolve('get_upbo', lp) */

static void impl_get_upbo(structlpsolve *lpsolve)
{
        Double  *dpr, *dpr0;

	if (lpsolve->lpsolvecaller.nrhs == 1+1) {
                int n, i;

                Check_nrhs(lpsolve, 1);
                n = get_Ncolumns(lpsolve->lp);
                dpr0 = dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, n, 1, 0);
                for (i = 1; i <= n; i++)
                  *(dpr++) = get_upbo(lpsolve->lp, i);
                SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr0, n, 1, 0, TRUE);
        }
        else {
	        Check_nrhs(lpsolve, 2);
                dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
		*dpr = get_upbo(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2));
                SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, 1, 0, TRUE);
        }
}


/* return = xxlpsolve('get_var_branch', lp, column) */
/* [var_branch] = xxlpsolve('get_var_branch', lp) */

static void impl_get_var_branch(structlpsolve *lpsolve)
{
        Long    *ipr, *ipr0;

	if (lpsolve->lpsolvecaller.nrhs == 1+1) {
                int i, n;

                Check_nrhs(lpsolve, 1);
                n = get_Ncolumns(lpsolve->lp);
                if (return_constants) {
                	char **types = NULL, *ptr, buf[bufsz];

                        types = (char **) callocmem(lpsolve, n, sizeof(*types));
                        for (i = 1; i <= n; i++) {
                                createconstant(lpsolve, get_var_branch(lpsolve->lp, i), consttype_branch, buf);
                                ptr = types[i - 1] = (char *) callocmem(lpsolve, strlen(buf) + 1, sizeof(**types));
                                strcpy(ptr, buf);
                        }
                        CreateString(&lpsolve->lpsolvecaller, types, n, 0);
                        for (i = 0; i < n; i++)
                                freemem(lpsolve, types[i]);
                        freemem(lpsolve, types);
                }
                else {
                        ipr0 = ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, n, 1, 0);
                        for (i = 1; i <= n; i++)
                                *(ipr++) = get_var_branch(lpsolve->lp, i);
                        SetLongMatrix(&lpsolve->lpsolvecaller, ipr0, n, 1, 0, TRUE);
                }
        }
        else {
        	Check_nrhs(lpsolve, 2);
                returnconstant(lpsolve, get_var_branch(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2)), consttype_branch);
        }
}


/* return = xxlpsolve('get_var_dualresult', lp, index) */

static void impl_get_var_dualresult(structlpsolve *lpsolve)
{
        Double  *dpr;

        Check_nrhs(lpsolve, 2);
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*dpr = get_var_dualresult(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2));
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_var_primalresult', lp, index) */

static void impl_get_var_primalresult(structlpsolve *lpsolve)
{
        Double  *dpr;

        Check_nrhs(lpsolve, 2);
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*dpr = get_var_primalresult(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2));
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('get_var_priority', lp, column) */
/* [var_priority] = xxlpsolve('get_var_priority', lp) */

static void impl_get_var_priority(structlpsolve *lpsolve)
{
        Long    *ipr, *ipr0;

	if (lpsolve->lpsolvecaller.nrhs == 1+1) {
                int n, i;

                Check_nrhs(lpsolve, 1);
                n = get_Ncolumns(lpsolve->lp);
                ipr0 = ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, n, 1, 0);
                for (i = 1; i <= n; i++)
                  *(ipr++) = get_var_priority(lpsolve->lp, i);
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr0, n, 1, 0, TRUE);
        }
        else {
        	Check_nrhs(lpsolve, 2);
                ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
		*ipr = get_var_priority(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2));
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
        }
}


/* [var, return] = xxlpsolve('get_variables', lp) */

static void impl_get_variables(structlpsolve *lpsolve)
{
        int n;
        int     result;
        Double  *dpr;
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        n = get_Ncolumns(lpsolve->lp);
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, n, 1, 0);
        result = get_variables(lpsolve->lp, dpr);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, n, 1, 0, TRUE);
        if (lpsolve->lpsolvecaller.nlhs > 1) {
                ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 1);
                *ipr = result;
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 1, TRUE);
        }
}


/* return = xxlpsolve('get_verbose', lp) */

static void impl_get_verbose(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 1);
        returnconstant(lpsolve, get_verbose(lpsolve->lp), consttype_verbose);
}


/* return = xxlpsolve('get_working_objective', lp) */

static void impl_get_working_objective(structlpsolve *lpsolve)
{
        Double  *dpr;

        Check_nrhs(lpsolve, 1);
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *dpr = get_working_objective(lpsolve->lp);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, 1, 0, TRUE);
}


/* [basisvector, return] = xxlpsolve('guess_basis', lp, [guessvector]) */

static void impl_guess_basis(structlpsolve *lpsolve)
{
        int i, n, m;
        REAL *guessvector;
        int *basisvector, *basisvector0;
        int     result;
        Long    *ipr, *ipr0;

        Check_nrhs(lpsolve, 2);

        n = get_Ncolumns(lpsolve->lp);
        m = get_Nrows(lpsolve->lp);
        guessvector = (REAL *) callocmem(lpsolve, 1 + n, sizeof(REAL));
        basisvector0 = basisvector = (int *) callocmem(lpsolve, 1 + n + m, sizeof(int));
	GetRealVector(&lpsolve->lpsolvecaller, 2, guessvector, 1, n, TRUE);
        result = guess_basis(lpsolve->lp, guessvector, basisvector);
        ipr0 = ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, n + m, 1, 0);
        for (i = 0; i < n + m; i++)
          *(ipr++) = *(++basisvector);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr0, n + m, 1, 0, TRUE);
        freemem(lpsolve, basisvector0);
        freemem(lpsolve, guessvector);
        if (lpsolve->lpsolvecaller.nlhs > 1) {
                ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 1);
                *ipr = result;
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 1, TRUE);
        }
}


/* return = xxlpsolve('has_BFP', lp) */

static void impl_has_BFP(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = has_BFP(lpsolve->lp);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('has_XLI', lp) */

static void impl_has_XLI(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = has_XLI(lpsolve->lp);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('is_add_rowmode', lp) */

static void impl_is_add_rowmode(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = is_add_rowmode(lpsolve->lp);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('is_anti_degen', lp, testmask) */

static void impl_is_anti_degen(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = is_anti_degen(lpsolve->lp, constant(lpsolve, 2, consttype_antidegen));
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('is_binary', lp, column) */
/* [binary] = xxlpsolve('is_binary', lp) */

static void impl_is_binary(structlpsolve *lpsolve)
{
        Long    *ipr, *ipr0;

        if (lpsolve->lpsolvecaller.nrhs == 1+1) {
                int n, i;

                Check_nrhs(lpsolve, 1);
                n = get_Ncolumns(lpsolve->lp);
                ipr0 = ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, n, 1, 0);
                for (i = 1; i <= n; i++)
        		*(ipr++) = is_binary(lpsolve->lp, i);
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr0, n, 1, 0, TRUE);
        }
        else {
        	Check_nrhs(lpsolve, 2);
                ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
		*ipr = is_binary(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2));
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
        }
}


/* return = xxlpsolve('is_break_at_first', lp) */

static void impl_is_break_at_first(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = is_break_at_first(lpsolve->lp);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('is_constr_type', lp, row, mask) */

static void impl_is_constr_type(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 3);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = is_constr_type(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2), constant(lpsolve, 3, consttype_constrainttype));
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('is_debug', lp) */

static void impl_is_debug(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = is_debug(lpsolve->lp);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('is_feasible', lp, [values] {, threshold}) */

static void impl_is_feasible(structlpsolve *lpsolve)
{
        int i, n;
	REAL	*vec, threshold;
        int     result;
        Long    *ipr;

        if (lpsolve->lpsolvecaller.nrhs == 2+1)
                n = 2;
        else
                n = 3;
        Check_nrhs(lpsolve, n);
        i = get_Nrows(lpsolve->lp) + get_Ncolumns(lpsolve->lp);
        vec = (REAL *) callocmem(lpsolve, 1 + i, sizeof(REAL));
	GetRealVector(&lpsolve->lpsolvecaller, 2, vec, 1, i, TRUE);
        if (n == 2)
                threshold = get_epsint(lpsolve->lp);
        else
                threshold = GetRealScalar(&lpsolve->lpsolvecaller, 3);
        result = is_feasible(lpsolve->lp, vec, threshold);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
	freemem(lpsolve, vec);
}


/* return = xxlpsolve('is_free', lp, column) */
/* return = xxlpsolve('is_unbounded', lp, column) */
/* [free] = xxlpsolve('is_free', lp) */
/* [free] = xxlpsolve('is_unbounded', lp) */

static void impl_is_free(structlpsolve *lpsolve)
{
        Long    *ipr, *ipr0;

        if (lpsolve->lpsolvecaller.nrhs == 1+1) {
                int n, i;

                Check_nrhs(lpsolve, 1);
                n = get_Ncolumns(lpsolve->lp);
                ipr0 = ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, n, 1, 0);
                for (i = 1; i <= n; i++)
        		*(ipr++) = is_unbounded(lpsolve->lp, i);
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr0, n, 1, 0, TRUE);
        }
        else {
        	Check_nrhs(lpsolve, 2);
                ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
		*ipr = is_unbounded(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2));
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
        }
}


/* return = xxlpsolve('is_infinite', lp, value) */

static void impl_is_infinite(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = is_infinite(lpsolve->lp, GetRealScalar(&lpsolve->lpsolvecaller, 2));
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('is_int', lp, column) */
/* [int] = xxlpsolve('is_int', lp) */

static void impl_is_int(structlpsolve *lpsolve)
{
        Long    *ipr, *ipr0;

        if (lpsolve->lpsolvecaller.nrhs == 1+1) {
                int n, i;

                Check_nrhs(lpsolve, 1);
                n = get_Ncolumns(lpsolve->lp);
                ipr0 = ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, n, 1, 0);
                for (i = 1; i <= n; i++)
        		*(ipr++) = is_int(lpsolve->lp, i);
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr0, n, 1, 0, TRUE);
        }
        else {
                Check_nrhs(lpsolve, 2);
                ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        	*ipr = is_int(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2));
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
        }
}


/* return = xxlpsolve('is_integerscaling', lp) */

static void impl_is_integerscaling(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = is_integerscaling(lpsolve->lp);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('is_maxim', lp) */

static void impl_is_maxim(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = is_maxim(lpsolve->lp);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('is_nativeBFP', lp) */

static void impl_is_nativeBFP(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = is_nativeBFP(lpsolve->lp);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('is_nativeXLI', lp) */

static void impl_is_nativeXLI(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = is_nativeXLI(lpsolve->lp);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('is_negative', lp, column) */
/* [negative] = xxlpsolve('is_negative', lp) */

static void impl_is_negative(structlpsolve *lpsolve)
{
        Long    *ipr, *ipr0;

        if (lpsolve->lpsolvecaller.nrhs == 1+1) {
                int n, i;

                Check_nrhs(lpsolve, 1);
                n = get_Ncolumns(lpsolve->lp);
                ipr0 = ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, n, 1, 0);
                for (i = 1; i <= n; i++)
        		*(ipr++) = is_negative(lpsolve->lp, i);
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr0, n, 1, 0, TRUE);
        }
        else {
	        Check_nrhs(lpsolve, 2);
                ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
		*ipr = is_negative(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2));
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
        }
}


/* return = xxlpsolve('is_piv_mode', lp, testmask) */

static void impl_is_piv_mode(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = is_piv_mode(lpsolve->lp, constant(lpsolve, 2, consttype_price));
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('is_piv_rule', lp, rule) */

static void impl_is_piv_rule(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = is_piv_rule(lpsolve->lp, constant(lpsolve, 2, consttype_pricer));
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('is_presolve', lp, testmask) */

static void impl_is_presolve(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = is_presolve(lpsolve->lp, constant(lpsolve, 2, consttype_presolve));
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('is_scalemode', lp, testmask) */

static void impl_is_scalemode(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = is_scalemode(lpsolve->lp, constant(lpsolve, 2, consttype_scale));
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('is_scaletype', lp, scaletype) */

static void impl_is_scaletype(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = is_scaletype(lpsolve->lp, constant(lpsolve, 2, consttype_scale));
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('is_semicont', lp, column) */
/* [semicont] = xxlpsolve('is_semicont', lp) */

static void impl_is_semicont(structlpsolve *lpsolve)
{
        Long    *ipr, *ipr0;

        if (lpsolve->lpsolvecaller.nrhs == 1+1) {
                int n, i;

                Check_nrhs(lpsolve, 1);
                n = get_Ncolumns(lpsolve->lp);
                ipr0 = ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, n, 1, 0);
                for (i = 1; i <= n; i++)
        		*(ipr++) = is_semicont(lpsolve->lp, i);
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr0, n, 1, 0, TRUE);
        }
        else {
        	Check_nrhs(lpsolve, 2);
                ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
		*ipr = is_semicont(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2));
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
        }
}


/* return = xxlpsolve('is_SOS_var', lp, column) */
/* [SOS_var] = xxlpsolve('is_SOS_var', lp) */

static void impl_is_SOS_var(structlpsolve *lpsolve)
{
        Long    *ipr, *ipr0;

        if (lpsolve->lpsolvecaller.nrhs == 1+1) {
                int n, i;

                Check_nrhs(lpsolve, 1);
                n = get_Ncolumns(lpsolve->lp);
                ipr0 = ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, n, 1, 0);
                for (i = 1; i <= n; i++)
        		*(ipr++) = is_SOS_var(lpsolve->lp, i);
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr0, n, 1, 0, TRUE);
        }
        else {
        	Check_nrhs(lpsolve, 2);
                ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
		*ipr = is_SOS_var(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2));
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
        }
}


/* return = xxlpsolve('is_trace', lp) */

static void impl_is_trace(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = is_nativeXLI(lpsolve->lp);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('is_use_names', lp, isrow) */

static void impl_is_use_names(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = (Long) is_use_names(lpsolve->lp, (MYBOOL) GetRealScalar(&lpsolve->lpsolvecaller, 2));
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* versionstring = xxlpsolve('lp_solve_version') */

static void impl_lp_solve_version(structlpsolve *lpsolve)
{
        int majorversion, minorversion, release, build;
        char buf[bufsz], *pbuf = buf;

        Check_nrhs(lpsolve, 0);
        lp_solve_version(&majorversion, &minorversion, &release, &build);
        sprintf(buf, "%d.%d.%d.%d", majorversion, minorversion, release, build);
	CreateString(&lpsolve->lpsolvecaller, &pbuf, 1, 0);
}


/* lp = xxlpsolve('make_lp', rows, columns) */

static void impl_make_lp(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = create_handle(lpsolve, make_lp((int) GetRealScalar(&lpsolve->lpsolvecaller, 1), (int) GetRealScalar(&lpsolve->lpsolvecaller, 2)), "make_lp failed");
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* lp = xxlpsolve('resize_lp', lp, rows, columns) */

static void impl_resize_lp(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 3);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = resize_lp(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2), (int) GetRealScalar(&lpsolve->lpsolvecaller, 3));
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* xxlpsolve('print_constraints', lp {, columns}) */

static void impl_print_constraints(structlpsolve *lpsolve)
{
        int n, columns;

        if (lpsolve->lpsolvecaller.nrhs == 1+1)
        	n = 1;
        else
        	n = 2;
        Check_nrhs(lpsolve, n);
        if (n == 1)
                columns = 1;
        else
                columns = (int) GetRealScalar(&lpsolve->lpsolvecaller, 2);
	print_constraints(lpsolve->lp, columns);
}


/* return = xxlpsolve('print_debugdump', lp, filename) */

static void impl_print_debugdump(structlpsolve *lpsolve)
{
        char filename[260];
        int     result;
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        GetString(&lpsolve->lpsolvecaller, NULL, 2, filename, sizeof(filename), TRUE);
        result = print_debugdump(lpsolve->lp, filename);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* xxlpsolve('print_duals', lp) */

static void impl_print_duals(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 1);
	print_duals(lpsolve->lp);
}


/* xxlpsolve('print_lp', lp) */

static void impl_print_lp(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 1);
	print_lp(lpsolve->lp);
}


/* xxlpsolve('print_objective', lp) */

static void impl_print_objective(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 1);
	print_objective(lpsolve->lp);
}


/* xxlpsolve('print_scales', lp) */

static void impl_print_scales(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 1);
	print_scales(lpsolve->lp);
}


/* xxlpsolve('print_solution', lp {, columns}) */

static void impl_print_solution(structlpsolve *lpsolve)
{
        int n, columns;

        if (lpsolve->lpsolvecaller.nrhs == 1+1)
        	n = 1;
        else
        	n = 2;
        Check_nrhs(lpsolve, n);
        if (n == 1)
                columns = 1;
        else
                columns = (int) GetRealScalar(&lpsolve->lpsolvecaller, 2);
	print_solution(lpsolve->lp, columns);
}


/* xxlpsolve('print_str', lp, str) */

static void impl_print_str(structlpsolve *lpsolve)
{
        char buf[bufsz];

        Check_nrhs(lpsolve, 2);
        GetString(&lpsolve->lpsolvecaller, NULL, 2, buf, bufsz, TRUE);
        print_str(lpsolve->lp, buf);
}


/* xxlpsolve('print_tableau', lp) */

static void impl_print_tableau(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 1);
	print_tableau(lpsolve->lp);
}


/* [handle_vec] = xxlpsolve('print_handle') */
/* print all used handles */

static void impl_print_handle(structlpsolve *lpsolve)
{
        int i, j, k, n;
        Long    *ipr, *ipr0;
        char size = FALSE;

	j = 0;
        for (i = 0; i <= lp_last; i++)
	  if (lp[i] != NULL)
	     j++;

        if (lpsolve->lpsolvecaller.nrhs == 1)
        	n = 0;
        else
        	n = 1;
        Check_nrhs(lpsolve, n);
        if (n == 0)
                size = FALSE;
        else
                size = (char) GetRealScalar(&lpsolve->lpsolvecaller, 1);

        if (size) {
                ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
                *ipr = j;
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
        }
        else {
                if (j)
                  k = 1;
                else
                  k = 0;
                ipr0 = ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, j, k, 0);
        	for (i = 0; i <= lp_last; i++)
        	  if (lp[i] != NULL)
                    *(ipr++) = i;
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr0, j, k, 0, TRUE);
        }
}


/* [handle_vec] = xxlpsolve('get_handle', 'name') */
/* get handle from model name */

static void impl_get_handle(structlpsolve *lpsolve)
{
        hashelem *hp;
        char buf[bufsz];
        Long    *ipr;

        Check_nrhs(lpsolve, 1);
        GetString(&lpsolve->lpsolvecaller, NULL, 1, buf, bufsz, TRUE);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        if (handlehash != NULL)
        	hp = findhash(buf, handlehash);
        else
                hp = NULL;
        if (hp == NULL)
                *ipr = -1;
        else
        	*ipr = hp->index;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* xxlpsolve('return_constants' [, asstring]) */

static void impl_return_constants(structlpsolve *lpsolve)
{
        Long    *ipr;

        if (lpsolve->lpsolvecaller.nrhs >= 1+1) {
                Check_nrhs(lpsolve, 1);
                return_constants = (char) GetRealScalar(&lpsolve->lpsolvecaller, 1);
        }
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = (int) return_constants;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* [return, info] = xxlpsolve('read_basis', lp, filename) */

static void impl_read_basis(structlpsolve *lpsolve)
{
        char filename[260];
        int     result;
        Long    *ipr;
#       define info filename

        Check_nrhs(lpsolve, 2);
        GetString(&lpsolve->lpsolvecaller, NULL, 2, filename, sizeof(filename), TRUE);
        result = read_basis(lpsolve->lp, filename, (lpsolve->lpsolvecaller.nlhs > 1) ? info : NULL);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
        if (lpsolve->lpsolvecaller.nlhs > 1) {
                char *ptr = info;
                CreateString(&lpsolve->lpsolvecaller, &ptr, 1, 1);
        }
#       undef info
}


/* lp = xxlpsolve('read_freeMPS', filename {, options}) */

static void impl_read_freeMPS(structlpsolve *lpsolve)
{
        int n, options;
        char filename[260];
        Long    *ipr;

        if (lpsolve->lpsolvecaller.nrhs == 1+1)
        	n = 1;
        else
        	n = 2;
        Check_nrhs(lpsolve, n);
        if (n >= 2)
        	options = constant(lpsolve, 2, consttype_verbose | consttype_MPS);
        else
                options = NORMAL;
	GetString(&lpsolve->lpsolvecaller, NULL, 1, filename, sizeof(filename), TRUE);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = create_handle(lpsolve, read_freeMPS(filename, options), "read_freeMPS can't read file.");
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* lp = xxlpsolve('read_lp_file', filename {, verbose {, lp_name}}) */
/* lp = xxlpsolve('read_lp', filename {, verbose {, lp_name}}) */
/* lp = xxlpsolve('read_LP', filename {, verbose {, lp_name}}) */

static void impl_read_LP(structlpsolve *lpsolve)
{
        int n, verbose;
        char filename[260], lp_name[50];
        Long    *ipr;

        if (lpsolve->lpsolvecaller.nrhs == 1+1)
        	n = 1;
        else if (lpsolve->lpsolvecaller.nrhs == 1+2)
        	n = 2;
        else
        	n = 3;
        Check_nrhs(lpsolve, n);
	GetString(&lpsolve->lpsolvecaller, NULL, 1, filename, sizeof(filename), TRUE);
        if (n >= 2)
                verbose = constant(lpsolve, 2, consttype_verbose);
        else
                verbose = NORMAL;
        if (n >= 3)
        	GetString(&lpsolve->lpsolvecaller, NULL, 3, lp_name, sizeof(lp_name), TRUE);
        else
                *lp_name = 0;
        lpsolve->lp = read_LP(filename, verbose, lp_name);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = create_handle(lpsolve, lpsolve->lp, "read_LP can't read file.");
        set_handlename(lpsolve->lp, lp_name, (int) *ipr);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* lp = xxlpsolve('read_mps', filename {, options}) */
/* lp = xxlpsolve('read_MPS', filename {, options}) */

static void impl_read_MPS(structlpsolve *lpsolve)
{
        int n, options;
        char filename[260], *name;
        Long    *ipr;

        if (lpsolve->lpsolvecaller.nrhs == 1+1)
        	n = 1;
        else
        	n = 2;
        Check_nrhs(lpsolve, n);
	GetString(&lpsolve->lpsolvecaller, NULL, 1, filename, sizeof(filename), TRUE);
        if (n >= 2)
        	options = constant(lpsolve, 2, consttype_verbose | consttype_MPS);
        else
                options = NORMAL;
        lpsolve->lp = read_MPS(filename, options);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = create_handle(lpsolve, lpsolve->lp, "read_MPS can't read file.");
        name = get_lp_name(lpsolve->lp);
        if (name != NULL)
                set_handlename(lpsolve->lp, name, (int) *ipr);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* lp = xxlpsolve('read_params', lp, filename {, options}) */

static void impl_read_params(structlpsolve *lpsolve)
{
        int n;
        char filename[260], options[50];
        Long    *ipr;

        if (lpsolve->lpsolvecaller.nrhs == 1+2)
        	n = 2;
        else
        	n = 3;
        Check_nrhs(lpsolve, n);
	GetString(&lpsolve->lpsolvecaller, NULL, 2, filename, sizeof(filename), TRUE);
        if (n >= 3)
        	GetString(&lpsolve->lpsolvecaller, NULL, 3, options, sizeof(options), TRUE);
        else
                *options = 0;
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = read_params(lpsolve->lp, filename, options);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* lp = xxlpsolve('read_XLI', xliname, modelname {, dataname {, options {, verbose}}} */

static void impl_read_XLI(structlpsolve *lpsolve)
{
        int n, verbose;
        char xliname[260], modelname[260], dataname[260], options[260];
        Long    *ipr;

        if (lpsolve->lpsolvecaller.nrhs == 1+2)
        	n = 2;
        else if (lpsolve->lpsolvecaller.nrhs == 1+3)
        	n = 3;
        else if (lpsolve->lpsolvecaller.nrhs == 1+4)
        	n = 4;
        else
        	n = 5;
        Check_nrhs(lpsolve, n);
	GetString(&lpsolve->lpsolvecaller, NULL, 1, xliname, sizeof(xliname), TRUE);
        GetString(&lpsolve->lpsolvecaller, NULL, 2, modelname, sizeof(modelname), TRUE);
        if (n >= 3)
        	GetString(&lpsolve->lpsolvecaller, NULL, 3, dataname, sizeof(dataname), TRUE);
        else
                *dataname = 0;
        if (n >= 4)
        	GetString(&lpsolve->lpsolvecaller, NULL, 4, options, sizeof(options), TRUE);
        else
                *options = 0;
        if (n >= 5)
        	verbose = constant(lpsolve, 5, consttype_verbose);
        else
                verbose = NORMAL;
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = create_handle(lpsolve, read_XLI(xliname, modelname, (*dataname) ? dataname : NULL, options, verbose), "read_XLI can't read file.");
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* xxlpsolve('reset_params', lp) */

static void impl_reset_params(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 1);
        reset_params(lpsolve->lp);
}


/* return = xxlpsolve('set_add_rowmode', lp, turnon) */

static void impl_set_add_rowmode(structlpsolve *lpsolve)
{
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = set_add_rowmode(lpsolve->lp, (MYBOOL) GetRealScalar(&lpsolve->lpsolvecaller, 2));
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* xxlpsolve('set_anti_degen', lp, anti_degen) */

static void impl_set_anti_degen(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_anti_degen(lpsolve->lp, constant(lpsolve, 2, consttype_antidegen));
}


/* return = xxlpsolve('set_basis', lp, [bascolumn], nonbasic) */

static void impl_set_basis(structlpsolve *lpsolve)
{
        int i, *bascolumn;
        MYBOOL nonbasic;
        int     result;
        Long    *ipr;

        Check_nrhs(lpsolve, 3);
        nonbasic = (MYBOOL) GetRealScalar(&lpsolve->lpsolvecaller, 3);
        i = get_Nrows(lpsolve->lp) + ((nonbasic) ? get_Ncolumns(lpsolve->lp) : 0);
        bascolumn = (int *) callocmem(lpsolve, 1 + i, sizeof(*bascolumn));
	GetIntVector(&lpsolve->lpsolvecaller, 2, bascolumn, 1, i, TRUE);
	result = set_basis(lpsolve->lp, bascolumn, nonbasic);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
	freemem(lpsolve, bascolumn);
}


/* xxlpsolve('set_basiscrash', lp, mode) */

static void impl_set_basiscrash(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_basiscrash(lpsolve->lp, constant(lpsolve, 2, consttype_crash));
}


/* xxlpsolve('set_basisvar', lp, basisPos, enteringCol) */

static void impl_set_basisvar(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 3);
	set_basisvar(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2), (int) GetRealScalar(&lpsolve->lpsolvecaller, 3));
}


/* xxlpsolve('set_bb_depthlimit', lp, bb_maxlevel) */

static void impl_set_bb_depthlimit(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_bb_depthlimit(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2));
}


/* xxlpsolve('set_bb_floorfirst', lp, bb_floorfirst) */

static void impl_set_bb_floorfirst(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_bb_floorfirst(lpsolve->lp, constant(lpsolve, 2, consttype_branch));
}


/* xxlpsolve('set_bb_rule', lp, bb_rule) */

static void impl_set_bb_rule(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_bb_rule(lpsolve->lp, constant(lpsolve, 2, consttype_node));
}


/* return = xxlpsolve('set_BFP', lp, filename) */

static void impl_set_BFP(structlpsolve *lpsolve)
{
        char filename[260];
        int     result;
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        GetString(&lpsolve->lpsolvecaller, NULL, 2, filename, sizeof(filename), TRUE);
        result = set_BFP(lpsolve->lp, filename);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('set_binary', lp, column, must_be_bin) */
/* return = xxlpsolve('set_binary', lp, [must_be_bin]) */

static void impl_set_binary(structlpsolve *lpsolve)
{
        int     result;
        Long    *ipr;

        if (lpsolve->lpsolvecaller.nrhs == 1+2) {
                int i, n, *vec;

                Check_nrhs(lpsolve, 2);
                n = get_Ncolumns(lpsolve->lp);
                vec = (int *) callocmem(lpsolve, n, sizeof(*vec));
        	GetIntVector(&lpsolve->lpsolvecaller, 2, vec, 0, n, TRUE);
                result = TRUE;
                for (i = 0; (i < n) && (result); i++)
                        result = set_binary(lpsolve->lp, i + 1, (MYBOOL) vec[i]);
                freemem(lpsolve, vec);
        }
        else {
                Check_nrhs(lpsolve, 3);
        	result = set_binary(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2), (MYBOOL) GetRealScalar(&lpsolve->lpsolvecaller, 3));
        }
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('set_bounds', lp, column, lower, upper) */
/* return = xxlpsolve('set_bounds', lp, [lower], [upper]) */

static void impl_set_bounds(structlpsolve *lpsolve)
{
        int     result;
        Long    *ipr;

        if (lpsolve->lpsolvecaller.nrhs == 1+3) {
                int i, n;
                REAL	*lower, *upper;

                Check_nrhs(lpsolve, 3);
                n = get_Ncolumns(lpsolve->lp);
                lower = (REAL *) callocmem(lpsolve, n, sizeof(REAL));
                upper = (REAL *) callocmem(lpsolve, n, sizeof(REAL));
        	GetRealVector(&lpsolve->lpsolvecaller, 2, lower, 0, n, TRUE);
                GetRealVector(&lpsolve->lpsolvecaller, 3, upper, 0, n, TRUE);
                result = TRUE;
                for (i = 0; (i < n) && (result); i++)
                        result = set_bounds(lpsolve->lp, i + 1, lower[i], upper[i]);
                freemem(lpsolve, upper);
                freemem(lpsolve, lower);
        }
        else {
	        Check_nrhs(lpsolve, 4);
		result = set_bounds(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2), GetRealScalar(&lpsolve->lpsolvecaller, 3), GetRealScalar(&lpsolve->lpsolvecaller, 4));
        }
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* xxlpsolve('set_bounds_tighter', lp, tighten) */

static void impl_set_bounds_tighter(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_bounds_tighter(lpsolve->lp, (MYBOOL) GetRealScalar(&lpsolve->lpsolvecaller, 2));
}


/* xxlpsolve('set_break_at_first', lp, break_at_first) */

static void impl_set_break_at_first(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_break_at_first(lpsolve->lp, (MYBOOL) GetRealScalar(&lpsolve->lpsolvecaller, 2));
}


/* xxlpsolve('set_break_at_value', lp, break_at_value) */

static void impl_set_break_at_value(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_break_at_value(lpsolve->lp, GetRealScalar(&lpsolve->lpsolvecaller, 2));
}


/* return = xxlpsolve('set_col_name', lp, column, name) */
/* return = xxlpsolve('set_col_name', lp, [names]) */

static void impl_set_col_name(structlpsolve *lpsolve)
{
        int     result;
        char buf[bufsz];
        Long    *ipr;

        if (lpsolve->lpsolvecaller.nrhs == 1+2) {
                int n, i;
                strArray pa;

                Check_nrhs(lpsolve, 2);
                n = get_Ncolumns(lpsolve->lp);
                pa = GetCellCharItems(&lpsolve->lpsolvecaller, 2, n, TRUE);
                result = TRUE;
                for (i = 0; (i < n) && (result); i++) {
                	GetCellString(&lpsolve->lpsolvecaller, pa, i, buf, bufsz);
                        result = set_col_name(lpsolve->lp, i + 1, buf);
                }
                FreeCellCharItems(pa, n);
        }
        else {
                Check_nrhs(lpsolve, 3);
                GetString(&lpsolve->lpsolvecaller, NULL, 3, buf, bufsz, TRUE);
                result = set_col_name(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2), buf);
        }
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* since always a matrix vector is provided to both set_column and set_columnex, always call the more
   performant sparse version of the two routines */
/* return = xxlpsolve('set_column', lp, col_no, [column]) */
/* return = xxlpsolve('set_columnex', lp, col_no, [column]) */

static void impl_set_column(structlpsolve *lpsolve)
{
        int m, count;
        int	*index;
	REAL	*vec;
        int     result;
        Long    *ipr;

        Check_nrhs(lpsolve, 3);
        m = get_Nrows(lpsolve->lp);
        vec = (REAL *) callocmem(lpsolve, 1 + m, sizeof(*vec));
        index = (int *) callocmem(lpsolve, 1 + m, sizeof(*index));
	count = GetRealSparseVector(&lpsolve->lpsolvecaller, 3, vec, index, 0, 1 + m, 0);
	result = set_columnex(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2), count, vec, index);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
        freemem(lpsolve, index);
        freemem(lpsolve, vec);
}


/* return = xxlpsolve('set_constr_type', lp, row, con_type) */
/* return = xxlpsolve('set_constr_type', lp, [con_type]) */

static void impl_set_constr_type(structlpsolve *lpsolve)
{
        int     result;
        Long    *ipr;

        if (lpsolve->lpsolvecaller.nrhs == 1+2) {
                int i, m, *vec;
                strArray pa;

                result = TRUE;

                Check_nrhs(lpsolve, 2);
                m = get_Nrows(lpsolve->lp);

                if ((pa = GetCellCharItems(&lpsolve->lpsolvecaller, 2, m, FALSE)) != NULL) {
                        char buf[bufsz];

                        for (i = 0; (i < m) && (result); i++) {
                        	GetCellString(&lpsolve->lpsolvecaller, pa, i, buf, bufsz);
                		result = set_constr_type(lpsolve->lp, i + 1, constantfromstr(lpsolve, buf, consttype_constrainttype));
                        }
                        FreeCellCharItems(pa, m);
                }
                else {
                        vec = (int *) callocmem(lpsolve, m, sizeof(*vec));
                	GetIntVector(&lpsolve->lpsolvecaller, 2, vec, 0, m, TRUE);
                        for (i = 0; (i < m) && (result); i++)
                                result = set_constr_type(lpsolve->lp, i + 1, vec[i]);
                        freemem(lpsolve, vec);
                }
        }
        else {
	        Check_nrhs(lpsolve, 3);
		result = set_constr_type(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2), constant(lpsolve, 3, consttype_constrainttype));
        }
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* xxlpsolve('set_debug', lp, debug) */

static void impl_set_debug(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_debug(lpsolve->lp, (MYBOOL) GetRealScalar(&lpsolve->lpsolvecaller, 2));
}


/* xxlpsolve('set_epsb', lp, epsb) */

static void impl_set_epsb(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_epsb(lpsolve->lp, GetRealScalar(&lpsolve->lpsolvecaller, 2));
}


/* xxlpsolve('set_epsd', lp, epsd) */

static void impl_set_epsd(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_epsd(lpsolve->lp, GetRealScalar(&lpsolve->lpsolvecaller, 2));
}


/* xxlpsolve('set_epsel', lp, epsel) */

static void impl_set_epsel(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_epsel(lpsolve->lp, GetRealScalar(&lpsolve->lpsolvecaller, 2));
}


/* xxlpsolve('set_epsint', lp, epsint) */

static void impl_set_epsint(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_epsint(lpsolve->lp, GetRealScalar(&lpsolve->lpsolvecaller, 2));
}


/* xxlpsolve('set_epslevel', lp, epslevel) */

static void impl_set_epslevel(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_epslevel(lpsolve->lp, constant(lpsolve, 2, consttype_eps));
}


/* xxlpsolve('set_epsperturb', lp, epsperturb) */

static void impl_set_epsperturb(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_epsperturb(lpsolve->lp, GetRealScalar(&lpsolve->lpsolvecaller, 2));
}


/* xxlpsolve('set_epspivot', lp, epspivot) */

static void impl_set_epspivot(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_epspivot(lpsolve->lp, GetRealScalar(&lpsolve->lpsolvecaller, 2));
}


/* return = xxlpsolve('set_free', lp, column) */
/* return = xxlpsolve('set_unbounded', lp, column) */

static void impl_set_free(structlpsolve *lpsolve)
{
        int     result;
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
	result = set_unbounded(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2));
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* xxlpsolve('set_improve', lp, improve) */

static void impl_set_improve(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_improve(lpsolve->lp, constant(lpsolve, 2, consttype_improve));
}


/* xxlpsolve('set_infinite', lp, infinite) */

static void impl_set_infinite(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_infinite(lpsolve->lp, GetRealScalar(&lpsolve->lpsolvecaller, 2));
}


/* return = xxlpsolve('set_int', lp, column, must_be_int) */
/* return = xxlpsolve('set_int', lp, [must_be_int]) */

static void impl_set_int(structlpsolve *lpsolve)
{
        int     result;
        Long    *ipr;

        if (lpsolve->lpsolvecaller.nrhs == 1+2) {
                int i, n, *vec;

                Check_nrhs(lpsolve, 2);
                n = get_Ncolumns(lpsolve->lp);
                vec = (int *) callocmem(lpsolve, n, sizeof(*vec));
        	GetIntVector(&lpsolve->lpsolvecaller, 2, vec, 0, n, TRUE);
                result = TRUE;
                for (i = 0; (i < n) && (result); i++)
                        result = set_int(lpsolve->lp, i + 1, (MYBOOL) vec[i]);
                freemem(lpsolve, vec);
        }
        else {
	        Check_nrhs(lpsolve, 3);
		result = set_int(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2), (MYBOOL) GetRealScalar(&lpsolve->lpsolvecaller, 3));
        }
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('set_lowbo', lp, column, value) */
/* return = xxlpsolve('set_lowbo', lp, [values]) */

static void impl_set_lowbo(structlpsolve *lpsolve)
{
        int     result;
        Long    *ipr;

        if (lpsolve->lpsolvecaller.nrhs == 1+2) {
                int i, n;
                REAL *vec;

                Check_nrhs(lpsolve, 2);
                n = get_Ncolumns(lpsolve->lp);
                vec = (REAL *) callocmem(lpsolve, n, sizeof(*vec));
        	GetRealVector(&lpsolve->lpsolvecaller, 2, vec, 0, n, TRUE);
                result = TRUE;
                for (i = 0; (i < n) && (result); i++)
                        result = set_lowbo(lpsolve->lp, i + 1, vec[i]);
                freemem(lpsolve, vec);
        }
        else {
	        Check_nrhs(lpsolve, 3);
		result = set_lowbo(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2), GetRealScalar(&lpsolve->lpsolvecaller, 3));
        }
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('set_lp_name', lp, name) */

static void impl_set_lp_name(structlpsolve *lpsolve)
{
        int     result;
        char buf[bufsz];
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        GetString(&lpsolve->lpsolvecaller, NULL, 2, buf, bufsz, TRUE);
        set_handlename(lpsolve->lp, buf, lpsolve->h);
        result = set_lp_name(lpsolve->lp, buf);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('set_mat', lp, [matrix]) */
/* return = xxlpsolve('set_mat', lp, row, column, value) */

static void impl_set_mat(structlpsolve *lpsolve)
{
        int     result;
        Long    *ipr;

        if (lpsolve->lpsolvecaller.nrhs == 1+2) {
                int m, n, j, *index, *index1, count;
                REAL *obj = NULL, *obj1 = NULL, *vec, *vec1, a;
                rMatrix mat;

                mat = GetpMatrix(&lpsolve->lpsolvecaller, 2);
		/* Called with a matrix argument */
		m = GetM(&lpsolve->lpsolvecaller, mat);
		n = GetN(&lpsolve->lpsolvecaller, mat);

/* Printf("get_Nrows = %d, get_Ncolumns = %d, m = %d, n = %d\n",
      get_Nrows(lpsolve->lp), get_Ncolumns(lpsolve->lp), m, n); */

		if ((get_Nrows(lpsolve->lp) != m) || (get_Ncolumns(lpsolve->lp) != n))
			ErrMsgTxt(&lpsolve->lpsolvecaller, "Invalid matrix dimension.");

                obj = obj1 = (REAL *) callocmem(lpsolve, 1 + n, sizeof(*obj));
                result = get_row(lpsolve->lp, 0, obj);
                vec = (REAL *) callocmem(lpsolve, 1 + m, sizeof(*vec));
                index = (int *) callocmem(lpsolve, 1 + m, sizeof(*index));
                for (j = 1; (j <= n) && (result); j++) {
                        vec1 = vec;
                        index1 = index;
                        count = 0;
                        if ((a = (*(++obj1))) != 0.0) {
                                *(vec1++) = a;
                                *(index1++) = 0;
                                count ++;
                        }
			count += GetRealSparseVector(&lpsolve->lpsolvecaller, 2, vec1, index1, 1, m, j);
                        result = set_columnex(lpsolve->lp, j, count, vec, index);
                }
                freemem(lpsolve, index);
                freemem(lpsolve, vec);
                freemem(lpsolve, obj);
                Check_nrhs(lpsolve, 2);
        }
        else { /* called with a single matrix element */
        	Check_nrhs(lpsolve, 4);
		result = set_mat(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2), (int) GetRealScalar(&lpsolve->lpsolvecaller, 3), GetRealScalar(&lpsolve->lpsolvecaller, 4));
        }
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* xxlpsolve('set_maxim', lp) */

static void impl_set_maxim(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 1);
        set_maxim(lpsolve->lp);
}


/* xxlpsolve('set_maxpivot', max_num_inv) */

static void impl_set_maxpivot(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_maxpivot(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2));
}


/* xxlpsolve('set_minim', lp) */

static void impl_set_minim(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 1);
        set_minim(lpsolve->lp);
}


/* xxlpsolve('set_mip_gap', lp, absolute, mip_gap) */

static void impl_set_mip_gap(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 3);
	set_mip_gap(lpsolve->lp, (MYBOOL) GetRealScalar(&lpsolve->lpsolvecaller, 2), GetRealScalar(&lpsolve->lpsolvecaller, 3));
}


/* xxlpsolve('set_negrange', negrange) */

static void impl_set_negrange(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_negrange(lpsolve->lp, GetRealScalar(&lpsolve->lpsolvecaller, 2));
}


/* return = xxlpsolve('set_obj', lp, column, value) */
/* return = xxlpsolve('set_obj', lp, [values]) */

static void impl_set_obj(structlpsolve *lpsolve)
{
        if (lpsolve->lpsolvecaller.nrhs == 1+2) {
                impl_set_obj_fn(lpsolve);
        }
        else {
                int     result;
                Long    *ipr;

	        Check_nrhs(lpsolve, 3);
		result = set_obj(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2), GetRealScalar(&lpsolve->lpsolvecaller, 3));
                ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	        *ipr = result;
                SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
        }
}


/* xxlpsolve('set_obj_bound', lp, obj_bound) */

static void impl_set_obj_bound(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_obj_bound(lpsolve->lp, GetRealScalar(&lpsolve->lpsolvecaller, 2));
}


/* since always a matrix vector is provided to both set_obj_fn and set_obj_fnex, always call the more
   performant sparse version of the two routines */
/* return = xxlpsolve('set_obj_fn', lp, [row]) */
/* return = xxlpsolve('set_obj_fnex', lp, [row]) */

static void impl_set_obj_fn(structlpsolve *lpsolve)
{
        int n, count;
        int	*index;
	REAL	*vec;
        int     result;
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        n = get_Ncolumns(lpsolve->lp);
        vec = (REAL *) callocmem(lpsolve, 1 + n, sizeof(*vec));
        index = (int *) callocmem(lpsolve, 1 + n, sizeof(*index));
	count = GetRealSparseVector(&lpsolve->lpsolvecaller, 2, vec, index, 1, n, 0);
	result = set_obj_fnex(lpsolve->lp, count, vec, index);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
        freemem(lpsolve, index);
        freemem(lpsolve, vec);
}


/* return = xxlpsolve('set_outputfile', lp, filename) */

static void impl_set_outputfile(structlpsolve *lpsolve)
{
        char filename[260];
        int     result;
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        GetString(&lpsolve->lpsolvecaller, NULL, 2, filename, sizeof(filename), TRUE);
        result = set_outputfile(lpsolve->lp, (*filename) ? filename : NULL);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* xxlpsolve('set_pivoting', lp, pivoting) */

static void impl_set_pivoting(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_pivoting(lpsolve->lp, constant(lpsolve, 2, consttype_pricer|consttype_price));
}


/* xxlpsolve('set_preferdual', lp, dodual) */

static void impl_set_preferdual(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_preferdual(lpsolve->lp, (MYBOOL) GetRealScalar(&lpsolve->lpsolvecaller, 2));
}


/* xxlpsolve('set_presolve', lp, do_presolve {, maxloops}) */

static void impl_set_presolve(structlpsolve *lpsolve)
{
        int maxloops;

        if (lpsolve->lpsolvecaller.nrhs == 1+2) {
                Check_nrhs(lpsolve, 2);
                maxloops = get_presolveloops(lpsolve->lp);
        }
        else {
                Check_nrhs(lpsolve, 3);
                maxloops = (int) GetRealScalar(&lpsolve->lpsolvecaller, 3);
        }
	set_presolve(lpsolve->lp, constant(lpsolve, 2, consttype_presolve), maxloops);
}


/* xxlpsolve('set_print_sol', lp, print_sol) */

static void impl_set_print_sol(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_print_sol(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2));
}


/* return = xxlpsolve('set_rh', lp, row, value) */
/* return = xxlpsolve('set_rh', lp, [values]) */

static void impl_set_rh(structlpsolve *lpsolve)
{
        int     result;
        Long    *ipr;

        if (lpsolve->lpsolvecaller.nrhs == 1+2) {
                int i, m;
                REAL *vec;

                Check_nrhs(lpsolve, 2);
                m = get_Nrows(lpsolve->lp);
                vec = (REAL *) callocmem(lpsolve, 1 + m, sizeof(*vec));
        	GetRealVector(&lpsolve->lpsolvecaller, 2, vec, 0, 1 + m, TRUE);
                result = TRUE;
                for (i = 0; (i <= m) && (result); i++)
                        result = set_rh(lpsolve->lp, i, vec[i]);
                freemem(lpsolve, vec);
        }
        else {
	        Check_nrhs(lpsolve, 3);
		result = set_rh(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2), GetRealScalar(&lpsolve->lpsolvecaller, 3));
        }
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('set_rh_range', lp, row, deltavalue) */
/* return = xxlpsolve('set_rh_range', lp, [deltavalues]) */

static void impl_set_rh_range(structlpsolve *lpsolve)
{
        int     result;
        Long    *ipr;

        if (lpsolve->lpsolvecaller.nrhs == 1+2) {
                int i, m;
                REAL *vec;

                Check_nrhs(lpsolve, 2);
                m = get_Nrows(lpsolve->lp);
                vec = (REAL *) callocmem(lpsolve, 1 + m, sizeof(*vec));
        	GetRealVector(&lpsolve->lpsolvecaller, 2, vec, 0, 1 + m, TRUE);
                result = TRUE;
                for (i = 0; (i < m) && (result); i++)
                        result = set_rh_range(lpsolve->lp, i + 1, vec[i]);
                freemem(lpsolve, vec);
        }
        else {
	        Check_nrhs(lpsolve, 3);
		result = set_rh_range(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2), GetRealScalar(&lpsolve->lpsolvecaller, 3));
        }
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* xxlpsolve('set_rh_vec', lp, [rh]) */

static void impl_set_rh_vec(structlpsolve *lpsolve)
{
        int m;
	REAL	*vec;

        Check_nrhs(lpsolve, 2);
        m = get_Nrows(lpsolve->lp);
	vec = (REAL *) callocmem(lpsolve, 1 + m, sizeof(REAL));
	GetRealVector(&lpsolve->lpsolvecaller, 2, vec, 1, m, TRUE);
	set_rh_vec(lpsolve->lp, vec);
        freemem(lpsolve, vec);
}


/* since always a matrix vector is provided to both set_row and set_rowex, always call the more
   performant sparse version of the two routines */
/* return = xxlpsolve('set_row', lp, row_no, [row]) */
/* return = xxlpsolve('set_rowex', lp, row_no, [row]) */

static void impl_set_row(structlpsolve *lpsolve)
{
        int n, count;
        int	*index;
	REAL	*vec;
        int     result;
        Long    *ipr;

        Check_nrhs(lpsolve, 3);
        n = get_Ncolumns(lpsolve->lp);
        vec = (REAL *) callocmem(lpsolve, 1 + n, sizeof(*vec));
        index = (int *) callocmem(lpsolve, 1 + n, sizeof(*index));
	count = GetRealSparseVector(&lpsolve->lpsolvecaller, 3, vec, index, 1, n, 0);
	result = set_rowex(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2), count, vec, index);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
        freemem(lpsolve, index);
        freemem(lpsolve, vec);
}


/* return = xxlpsolve('set_row_name', lp, row, name) */
/* return = xxlpsolve('set_row_name', lp, [names]) */

static void impl_set_row_name(structlpsolve *lpsolve)
{
        int     result;
        char buf[bufsz];
        Long    *ipr;

        if (lpsolve->lpsolvecaller.nrhs == 1+2) {
                int m, i;
                strArray pa;

                Check_nrhs(lpsolve, 2);
                m = get_Nrows(lpsolve->lp);
                pa = GetCellCharItems(&lpsolve->lpsolvecaller, 2, m, TRUE);
                result = TRUE;
                for (i = 0; (i < m) && (result); i++) {
                	GetCellString(&lpsolve->lpsolvecaller, pa, i, buf, bufsz);
                        result = set_row_name(lpsolve->lp, i + 1, buf);
                }
                FreeCellCharItems(pa, m);
        }
        else {
	        Check_nrhs(lpsolve, 3);
	        GetString(&lpsolve->lpsolvecaller, NULL, 3, buf, bufsz, TRUE);
	        result = set_row_name(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2), buf);
        }
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* xxlpsolve('set_scalelimit', lp, scalelimit) */

static void impl_set_scalelimit(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_scalelimit(lpsolve->lp, GetRealScalar(&lpsolve->lpsolvecaller, 2));
}


/* xxlpsolve('set_scaling', lp, scalemode) */

static void impl_set_scaling(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_scaling(lpsolve->lp, constant(lpsolve, 2, consttype_scale));
}


/* return = xxlpsolve('set_semicont', lp, column, must_be_sc) */
/* return = xxlpsolve('set_semicont', lp, [must_be_sc]) */

static void impl_set_semicont(structlpsolve *lpsolve)
{
        int     result;
        Long    *ipr;

        if (lpsolve->lpsolvecaller.nrhs == 1+2) {
                int i, n, *vec;

                Check_nrhs(lpsolve, 2);
                n = get_Ncolumns(lpsolve->lp);
                vec = (int *) callocmem(lpsolve, n, sizeof(*vec));
        	GetIntVector(&lpsolve->lpsolvecaller, 2, vec, 0, n, TRUE);
                result = TRUE;
                for (i = 0; (i < n) && (result); i++)
                        result = set_semicont(lpsolve->lp, i + 1, (MYBOOL) vec[i]);
                freemem(lpsolve, vec);
        }
        else {
	        Check_nrhs(lpsolve, 3);
		result = set_semicont(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2), (MYBOOL) GetRealScalar(&lpsolve->lpsolvecaller, 3));
        }
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* xxlpsolve('set_sense', lp, maximize) */

static void impl_set_sense(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_sense(lpsolve->lp, (MYBOOL) GetRealScalar(&lpsolve->lpsolvecaller, 2));
}


/* xxlpsolve('set_simplextype', lp, simplextype) */

static void impl_set_simplextype(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_simplextype(lpsolve->lp, constant(lpsolve, 2, consttype_simplex));
}


/* xxlpsolve('set_solutionlimit', lp, simplextype) */

static void impl_set_solutionlimit(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_solutionlimit(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2));
}


/* xxlpsolve('set_timeout', lp, sectimeout) */

static void impl_set_timeout(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_timeout(lpsolve->lp, (long) GetRealScalar(&lpsolve->lpsolvecaller, 2));
}


/* xxlpsolve('set_trace', lp, trace) */

static void impl_set_trace(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_trace(lpsolve->lp, (MYBOOL) GetRealScalar(&lpsolve->lpsolvecaller, 2));
}


/* return = xxlpsolve('set_upbo', lp, column, value) */
/* return = xxlpsolve('set_upbo', lp, [values]) */

static void impl_set_upbo(structlpsolve *lpsolve)
{
        int     result;
        Long    *ipr;

        if (lpsolve->lpsolvecaller.nrhs == 1+2) {
                int i, n;
                REAL *vec;

                Check_nrhs(lpsolve, 2);
                n = get_Ncolumns(lpsolve->lp);
                vec = (REAL *) callocmem(lpsolve, n, sizeof(*vec));
        	GetRealVector(&lpsolve->lpsolvecaller, 2, vec, 0, n, TRUE);
                result = TRUE;
                for (i = 0; (i < n) && (result); i++)
                        result = set_upbo(lpsolve->lp, i + 1, vec[i]);
                freemem(lpsolve, vec);
        }
        else {
	        Check_nrhs(lpsolve, 3);
		result = set_upbo(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2), GetRealScalar(&lpsolve->lpsolvecaller, 3));
        }
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* xxlpsolve('set_use_names', lp, isrow, use_names) */

static void impl_set_use_names(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 3);
        set_use_names(lpsolve->lp, (MYBOOL) GetRealScalar(&lpsolve->lpsolvecaller, 2), (MYBOOL) GetRealScalar(&lpsolve->lpsolvecaller, 3));
}


/* return = xxlpsolve('set_var_branch', lp, column, branch_mode) */
/* return = xxlpsolve('set_var_branch', lp, [branch_mode]) */

static void impl_set_var_branch(structlpsolve *lpsolve)
{
        int     result;
        Long    *ipr;

        if (lpsolve->lpsolvecaller.nrhs == 1+2) {
                int i, n, *vec;
                strArray pa;

                result = TRUE;

                Check_nrhs(lpsolve, 2);
                n = get_Ncolumns(lpsolve->lp);

                if ((pa = GetCellCharItems(&lpsolve->lpsolvecaller, 2, n, FALSE)) != NULL) {
                        char buf[bufsz];

                        for (i = 0; (i < n) && (result); i++) {
                        	GetCellString(&lpsolve->lpsolvecaller, pa, i, buf, bufsz);
                		result = set_var_branch(lpsolve->lp, i + 1, constantfromstr(lpsolve, buf, consttype_branch));
                        }
                        FreeCellCharItems(pa, n);
                }
                else {
                        vec = (int *) callocmem(lpsolve, n, sizeof(*vec));
                	GetIntVector(&lpsolve->lpsolvecaller, 2, vec, 0, n, TRUE);
                        for (i = 0; (i < n) && (result); i++)
                                result = set_var_branch(lpsolve->lp, i + 1, vec[i]);
                        freemem(lpsolve, vec);
                }
        }
        else {
	        Check_nrhs(lpsolve, 3);
		result = set_var_branch(lpsolve->lp, (int) GetRealScalar(&lpsolve->lpsolvecaller, 2), constant(lpsolve, 3, consttype_branch));
        }
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('set_var_weights', lp, [weights]) */

static void impl_set_var_weights(structlpsolve *lpsolve)
{
        int n;
	REAL	*vec;
        int     result;
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        n = get_Ncolumns(lpsolve->lp);
	vec = (REAL *) callocmem(lpsolve, n, sizeof(REAL));
	GetRealVector(&lpsolve->lpsolvecaller, 2, vec, 0, n, TRUE);
	result = set_var_weights(lpsolve->lp, vec);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
        freemem(lpsolve, vec);
}


/* xxlpsolve('set_verbose', lp, verbose) */

static void impl_set_verbose(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 2);
	set_verbose(lpsolve->lp, constant(lpsolve, 2, consttype_verbose));
}


/* return = xxlpsolve('set_XLI', lp, filename) */

static void impl_set_XLI(structlpsolve *lpsolve)
{
        char filename[260];
        int     result;
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        GetString(&lpsolve->lpsolvecaller, NULL, 2, filename, sizeof(filename), TRUE);
        result = set_XLI(lpsolve->lp, filename);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* result = xxlpsolve('solve', lp) */

static void impl_solve(structlpsolve *lpsolve)
{
        /* int m, n, i ; */
        int     result;

        Check_nrhs(lpsolve, 1);
	result = solve(lpsolve->lp);
        returnconstant(lpsolve, result, consttype_solve);
        switch (result) {
        case OPTIMAL:
        case SUBOPTIMAL:
        case PROCBREAK:
        case FEASFOUND:
        case PRESOLVED:
/*
          if (get_verbose(lpsolve->lp) >= DETAILED) {
	    Printf("Branch & Bound depth: %d%s", get_max_level(lpsolve->lp), NEWLINE);
	    Printf("Nodes processed: %.0f%s", (double) get_total_nodes(lpsolve->lp), NEWLINE);
	    Printf("Simplex pivots: %.0f%s", (double) get_total_iter(lpsolve->lp), NEWLINE);
            Printf("Number of equal solutions: %d%s%s", get_solutioncount(lpsolve->lp), NEWLINE, NEWLINE);
          }
          if (get_verbose(lpsolve->lp) >= NORMAL) {
            Printf("Value of objective function: %f%s%s", get_objective(lpsolve->lp), NEWLINE, NEWLINE);
            Printf("Actual values of the variables:%s", NEWLINE);
            m = get_Nrows(lpsolve->lp);
            n = get_Ncolumns(lpsolve->lp);
            for (i = 1; i <= n; i++)
              Printf("%s: %f%s", get_col_name(lpsolve->lp, i), get_var_primalresult(lpsolve->lp, m + i), NEWLINE);
          }
*/
          break;
        case NOMEMORY:
          if (get_verbose(lpsolve->lp) >= NORMAL)
          	Printf("Out of memory%s", NEWLINE);
          break;
        case INFEASIBLE:
          if (get_verbose(lpsolve->lp) >= NORMAL)
	  	Printf("This problem is infeasible%s", NEWLINE);
          break;
        case UNBOUNDED:
          if (get_verbose(lpsolve->lp) >= NORMAL)
          	Printf("This problem is unbounded%s", NEWLINE);
          break;
        case PROCFAIL:
          if (get_verbose(lpsolve->lp) >= NORMAL)
          	Printf("The B&B routine failed%s", NEWLINE);
          break;
        case TIMEOUT:
          if (get_verbose(lpsolve->lp) >= NORMAL)
          	Printf("Timeout%s", NEWLINE);
          break;
        case USERABORT:
          if (get_verbose(lpsolve->lp) >= NORMAL)
          	Printf("User aborted%s", NEWLINE);
          break;
        case DEGENERATE:
          if (get_verbose(lpsolve->lp) >= NORMAL)
          	Printf("This problem is degenerative%s", NEWLINE);
          break;
        case NUMFAILURE:
          if (get_verbose(lpsolve->lp) >= NORMAL)
          	Printf("Numerical failure encountered%s", NEWLINE);
          break;
        case NOFEASFOUND:
          if (get_verbose(lpsolve->lp) >= NORMAL)
          	Printf("No feasible branch and bound solution found%s", NEWLINE);
          break;
        default:
          if (get_verbose(lpsolve->lp) >= NORMAL)
          	Printf("lp_solve failed%s", NEWLINE);
          break;
        }
}


/* return = xxlpsolve('time_elapsed', lp) */

static void impl_time_elapsed(structlpsolve *lpsolve)
{
        Double  *dpr;

        Check_nrhs(lpsolve, 1);
        dpr = CreateDoubleMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *dpr = time_elapsed(lpsolve->lp);
        SetDoubleMatrix(&lpsolve->lpsolvecaller, dpr, 1, 1, 0, TRUE);
}


/* xxlpsolve('unscale', lp) */

static void impl_unscale(structlpsolve *lpsolve)
{
        Check_nrhs(lpsolve, 1);
        unscale(lpsolve->lp);
}


/* return = xxlpsolve('write_basis', lp, filename) */

static void impl_write_basis(structlpsolve *lpsolve)
{
        char filename[260];
        int     result;
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        GetString(&lpsolve->lpsolvecaller, NULL, 2, filename, sizeof(filename), TRUE);
        result = write_basis(lpsolve->lp, filename);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('write_freemps', lp, filename) */
/* return = xxlpsolve('write_freeMPS', lp, filename) */

static void impl_write_freemps(structlpsolve *lpsolve)
{
        char filename[260];
        int     result;
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        GetString(&lpsolve->lpsolvecaller, NULL, 2, filename, sizeof(filename), TRUE);
        result = write_freemps(lpsolve->lp, filename);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('write_lp', lp, filename) */
/* return = xxlpsolve('write_LP', lp, filename) */

static void impl_write_lp(structlpsolve *lpsolve)
{
        char filename[260];
        int     result;
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        GetString(&lpsolve->lpsolvecaller, NULL, 2, filename, sizeof(filename), TRUE);
        result = write_lp(lpsolve->lp, filename);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('write_mps', lp, filename) */
/* return = xxlpsolve('write_MPS', lp, filename) */

static void impl_write_mps(structlpsolve *lpsolve)
{
        char filename[260];
        int     result;
        Long    *ipr;

        Check_nrhs(lpsolve, 2);
        GetString(&lpsolve->lpsolvecaller, NULL, 2, filename, sizeof(filename), TRUE);
        result = write_mps(lpsolve->lp, filename);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* lp = xxlpsolve('write_params', lp, filename {, options}) */

static void impl_write_params(structlpsolve *lpsolve)
{
        int n;
        char filename[260], options[50];
        Long    *ipr;

        if (lpsolve->lpsolvecaller.nrhs == 1+2)
        	n = 2;
        else
        	n = 3;
        Check_nrhs(lpsolve, n);
	GetString(&lpsolve->lpsolvecaller, NULL, 2, filename, sizeof(filename), TRUE);
        if (n >= 3)
        	GetString(&lpsolve->lpsolvecaller, NULL, 3, options, sizeof(options), TRUE);
        else
                *options = 0;
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
        *ipr = write_params(lpsolve->lp, filename, options);
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


/* return = xxlpsolve('write_XLI', lp, filename {, options {, results}}) */

static void impl_write_XLI(structlpsolve *lpsolve)
{
        char filename[260], options[50];
        int n;
        MYBOOL results;
        int     result;
        Long    *ipr;

        if (lpsolve->lpsolvecaller.nrhs == 1+2)
                n = 2;
        else if (lpsolve->lpsolvecaller.nrhs == 1+3)
                n = 3;
        else
                n = 4;

        Check_nrhs(lpsolve, n);
        GetString(&lpsolve->lpsolvecaller, NULL, 2, filename, sizeof(filename), TRUE);
        if (n >= 3)
        	GetString(&lpsolve->lpsolvecaller, NULL, 3, options, sizeof(options), TRUE);
        else
                *options = 0;
        if (n >= 4)
                results = (MYBOOL) GetRealScalar(&lpsolve->lpsolvecaller, 4);
        else
                results = FALSE;
        result = write_XLI(lpsolve->lp, filename, options, results);
        ipr = CreateLongMatrix(&lpsolve->lpsolvecaller, 1, 1, 0);
	*ipr = result;
        SetLongMatrix(&lpsolve->lpsolvecaller, ipr, 1, 1, 0, TRUE);
}


static struct {
        char *cmd;
        impl_routine *routine;
        int needshandle;
} routines[] = {
  { "add_column", impl_add_column, TRUE },
  { "add_columnex", impl_add_column, TRUE },
  { "add_constraint", impl_add_constraint, TRUE },
  { "add_constraintex", impl_add_constraint, TRUE },
  { "add_SOS", impl_add_SOS, TRUE },
  { "column_in_lp", impl_column_in_lp, TRUE },
  { "copy_lp", impl_copy_lp, TRUE },
  { "default_basis", impl_default_basis, TRUE },
  { "del_column", impl_del_column, TRUE },
  { "del_constraint", impl_del_constraint, TRUE },
  { "delete_lp", impl_delete_lp, TRUE },
  { "dualize_lp", impl_dualize_lp, TRUE },
  { "free_lp", impl_delete_lp, TRUE },
  { "get_anti_degen", impl_get_anti_degen, TRUE },
  { "get_basis", impl_get_basis, TRUE },
  { "get_basiscrash", impl_get_basiscrash, TRUE },
  { "get_bb_depthlimit", impl_get_bb_depthlimit, TRUE },
  { "get_bb_floorfirst", impl_get_bb_floorfirst, TRUE },
  { "get_bb_rule", impl_get_bb_rule, TRUE },
  { "get_bounds_tighter", impl_get_bounds_tighter, TRUE },
  { "get_break_at_value", impl_get_break_at_value, TRUE },
  { "get_col_name", impl_get_col_name, TRUE },
  { "get_column", impl_get_column, TRUE },
  { "get_columnex", impl_get_column, TRUE },
  { "get_constr_type", impl_get_constr_type, TRUE },
  { "get_constr_value", impl_get_constr_value, TRUE },
  { "get_constraints", impl_get_constraints, TRUE },
  { "get_dual_solution", impl_get_dual_solution, TRUE },
  { "get_epsb", impl_get_epsb, TRUE },
  { "get_epsd", impl_get_epsd, TRUE },
  { "get_epsel", impl_get_epsel, TRUE },
  { "get_epsint", impl_get_epsint, TRUE },
  { "get_epsperturb", impl_get_epsperturb, TRUE },
  { "get_epspivot", impl_get_epspivot, TRUE },
  { "get_improve", impl_get_improve, TRUE },
  { "get_infinite", impl_get_infinite, TRUE },
  { "get_lowbo", impl_get_lowbo, TRUE },
  { "get_lp_index", impl_get_lp_index, TRUE },
  { "get_lp_name", impl_get_lp_name, TRUE },
  { "get_mat", impl_get_mat, TRUE },
  { "get_max_level", impl_get_max_level, TRUE },
  { "get_maxpivot", impl_get_maxpivot, TRUE },
  { "get_mip_gap", impl_get_mip_gap, TRUE },
  { "get_nameindex", impl_get_nameindex, TRUE },
  { "get_Ncolumns", impl_get_Ncolumns, TRUE },
  { "get_negrange", impl_get_negrange, TRUE },
  { "get_nonzeros", impl_get_nonzeros, TRUE },
  { "get_Norig_columns", impl_get_Norig_columns, TRUE },
  { "get_Norig_rows", impl_get_Norig_rows, TRUE },
  { "get_Nrows", impl_get_Nrows, TRUE },
  { "get_obj_bound", impl_get_obj_bound, TRUE },
  { "get_objective", impl_get_objective, TRUE },
  { "get_orig_index", impl_get_orig_index, TRUE },
  { "get_origcol_name", impl_get_origcol_name, TRUE },
  { "get_origrow_name", impl_get_origrow_name, TRUE },
  { "get_pivoting", impl_get_pivoting, TRUE },
  { "get_presolve", impl_get_presolve, TRUE },
  { "get_presolveloops", impl_get_presolveloops, TRUE },
  { "get_primal_solution", impl_get_primal_solution, TRUE },
  { "get_print_sol", impl_get_print_sol, TRUE },
  { "get_rh", impl_get_rh, TRUE },
  { "get_rh_range", impl_get_rh_range, TRUE },
  { "get_row", impl_get_row, TRUE },
  { "get_rowex", impl_get_row, TRUE },
  { "get_row_name", impl_get_row_name, TRUE },
  { "get_scalelimit", impl_get_scalelimit, TRUE },
  { "get_scaling", impl_get_scaling, TRUE },
  { "get_sensitivity_obj", impl_get_sensitivity_objex, TRUE },
  { "get_sensitivity_objex", impl_get_sensitivity_objex, TRUE },
  { "get_sensitivity_rhs", impl_get_sensitivity_rhsex, TRUE },
  { "get_sensitivity_rhsex", impl_get_sensitivity_rhsex, TRUE },
  { "get_simplextype", impl_get_simplextype, TRUE },
  { "get_solution", impl_get_solution, TRUE },
  { "get_solutioncount", impl_get_solutioncount, TRUE },
  { "get_solutionlimit", impl_get_solutionlimit, TRUE },
  { "get_status", impl_get_status, TRUE },
  { "get_statustext", impl_get_statustext, TRUE },
  { "get_timeout", impl_get_timeout, TRUE },
  { "get_total_iter", impl_get_total_iter, TRUE },
  { "get_total_nodes", impl_get_total_nodes, TRUE },
  { "get_upbo", impl_get_upbo, TRUE },
  { "get_var_branch", impl_get_var_branch, TRUE },
  { "get_var_dualresult", impl_get_var_dualresult, TRUE },
  { "get_var_primalresult", impl_get_var_primalresult, TRUE },
  { "get_var_priority", impl_get_var_priority, TRUE },
  { "get_variables", impl_get_variables, TRUE },
  { "get_verbose", impl_get_verbose, TRUE },
  { "get_working_objective", impl_get_working_objective, TRUE },
  { "guess_basis", impl_guess_basis, TRUE },
  { "has_BFP", impl_has_BFP, TRUE },
  { "has_XLI", impl_has_XLI, TRUE },
  { "is_add_rowmode", impl_is_add_rowmode, TRUE },
  { "is_anti_degen", impl_is_anti_degen, TRUE },
  { "is_binary", impl_is_binary, TRUE },
  { "is_break_at_first", impl_is_break_at_first, TRUE },
  { "is_constr_type", impl_is_constr_type, TRUE },
  { "is_debug", impl_is_debug, TRUE },
  { "is_feasible", impl_is_feasible, TRUE },
  { "is_free", impl_is_free, TRUE },
  { "is_infinite", impl_is_infinite, TRUE },
  { "is_int", impl_is_int, TRUE },
  { "is_integerscaling", impl_is_integerscaling, TRUE },
  { "is_maxim", impl_is_maxim, TRUE },
  { "is_nativeBFP", impl_is_nativeBFP, TRUE },
  { "is_nativeXLI", impl_is_nativeXLI, TRUE },
  { "is_negative", impl_is_negative, TRUE },
  { "is_piv_mode", impl_is_piv_mode, TRUE },
  { "is_piv_rule", impl_is_piv_rule, TRUE },
  { "is_presolve", impl_is_presolve, TRUE },
  { "is_scalemode", impl_is_scalemode, TRUE },
  { "is_scaletype", impl_is_scaletype, TRUE },
  { "is_semicont", impl_is_semicont, TRUE },
  { "is_SOS_var", impl_is_SOS_var, TRUE },
  { "is_trace", impl_is_trace, TRUE },
  { "is_unbounded", impl_is_free, TRUE },
  { "is_use_names", impl_is_use_names, TRUE },
  { "lp_solve_version", impl_lp_solve_version, FALSE },
  { "make_lp", impl_make_lp, FALSE },
  { "print_constraints", impl_print_constraints, TRUE },
  { "print_debugdump", impl_print_debugdump, TRUE },
  { "print_duals", impl_print_duals, TRUE },
  { "print_lp", impl_print_lp, TRUE },
  { "print_objective", impl_print_objective, TRUE },
  { "print_scales", impl_print_scales, TRUE },
  { "print_solution", impl_print_solution, TRUE },
  { "print_str", impl_print_str, TRUE },
  { "print_tableau", impl_print_tableau, TRUE },
  { "read_basis", impl_read_basis, TRUE },
  { "read_freemps", impl_read_freeMPS, FALSE },
  { "read_freeMPS", impl_read_freeMPS, FALSE },
  { "read_lp", impl_read_LP, FALSE },
  { "read_LP", impl_read_LP, FALSE },
  { "read_mps", impl_read_MPS, FALSE },
  { "read_MPS", impl_read_MPS, FALSE },
  { "read_params", impl_read_params, TRUE },
  { "read_XLI", impl_read_XLI, FALSE },
  { "reset_params", impl_reset_params, TRUE },
  { "resize_lp", impl_resize_lp, TRUE },
  { "set_add_rowmode", impl_set_add_rowmode, TRUE },
  { "set_anti_degen", impl_set_anti_degen, TRUE },
  { "set_basis", impl_set_basis, TRUE },
  { "set_basiscrash", impl_set_basiscrash, TRUE },
  { "set_basisvar", impl_set_basisvar, TRUE },
  { "set_bb_depthlimit", impl_set_bb_depthlimit, TRUE },
  { "set_bb_floorfirst", impl_set_bb_floorfirst, TRUE },
  { "set_bb_rule", impl_set_bb_rule, TRUE },
  { "set_BFP", impl_set_BFP, TRUE },
  { "set_binary", impl_set_binary, TRUE },
  { "set_bounds", impl_set_bounds, TRUE },
  { "set_bounds_tighter", impl_set_bounds_tighter, TRUE },
  { "set_break_at_first", impl_set_break_at_first, TRUE },
  { "set_break_at_value", impl_set_break_at_value, TRUE },
  { "set_col_name", impl_set_col_name, TRUE },
  { "set_column", impl_set_column, TRUE },
  { "set_columnex", impl_set_column, TRUE },
  { "set_constr_type", impl_set_constr_type, TRUE },
  { "set_debug", impl_set_debug, TRUE },
  { "set_epsb", impl_set_epsb, TRUE },
  { "set_epsd", impl_set_epsd, TRUE },
  { "set_epsel", impl_set_epsel, TRUE },
  { "set_epsint", impl_set_epsint, TRUE },
  { "set_epslevel", impl_set_epslevel, TRUE },
  { "set_epsperturb", impl_set_epsperturb, TRUE },
  { "set_epspivot", impl_set_epspivot, TRUE },
  { "set_free", impl_set_free, TRUE },
  { "set_improve", impl_set_improve, TRUE },
  { "set_infinite", impl_set_infinite, TRUE },
  { "set_int", impl_set_int, TRUE },
  { "set_lowbo", impl_set_lowbo, TRUE },
  { "set_lp_name", impl_set_lp_name, TRUE },
  { "set_mat", impl_set_mat, TRUE },
  { "set_maxim", impl_set_maxim, TRUE },
  { "set_maxpivot", impl_set_maxpivot, TRUE },
  { "set_minim", impl_set_minim, TRUE },
  { "set_mip_gap", impl_set_mip_gap, TRUE },
  { "set_negrange", impl_set_negrange, TRUE },
  { "set_obj", impl_set_obj, TRUE },
  { "set_obj_bound", impl_set_obj_bound, TRUE },
  { "set_obj_fn", impl_set_obj_fn, TRUE },
  { "set_obj_fnex", impl_set_obj_fn, TRUE },
  { "set_outputfile", impl_set_outputfile, TRUE },
  { "set_pivoting", impl_set_pivoting, TRUE },
  { "set_preferdual", impl_set_preferdual, TRUE },
  { "set_presolve", impl_set_presolve, TRUE },
  { "set_print_sol", impl_set_print_sol, TRUE },
  { "set_rh", impl_set_rh, TRUE },
  { "set_rh_range", impl_set_rh_range, TRUE },
  { "set_rh_vec", impl_set_rh_vec, TRUE },
  { "set_row", impl_set_row, TRUE },
  { "set_rowex", impl_set_row, TRUE },
  { "set_row_name", impl_set_row_name, TRUE },
  { "set_scalelimit", impl_set_scalelimit, TRUE },
  { "set_scaling", impl_set_scaling, TRUE },
  { "set_semicont", impl_set_semicont, TRUE },
  { "set_sense", impl_set_sense, TRUE },
  { "set_simplextype", impl_set_simplextype, TRUE },
  { "set_solutionlimit", impl_set_solutionlimit, TRUE },
  { "set_timeout", impl_set_timeout, TRUE },
  { "set_trace", impl_set_trace, TRUE },
  { "set_unbounded", impl_set_free, TRUE },
  { "set_upbo", impl_set_upbo, TRUE },
  { "set_use_names", impl_set_use_names, TRUE },
  { "set_var_branch", impl_set_var_branch, TRUE },
  { "set_var_weights", impl_set_var_weights, TRUE },
  { "set_verbose", impl_set_verbose, TRUE },
  { "set_XLI", impl_set_XLI, TRUE },
  { "solve", impl_solve, TRUE },
  { "time_elapsed", impl_time_elapsed, TRUE },
  { "unscale", impl_unscale, TRUE },
  { "write_basis", impl_write_basis, TRUE },
  { "write_freemps", impl_write_freemps, TRUE },
  { "write_freeMPS", impl_write_freemps, TRUE },
  { "write_lp", impl_write_lp, TRUE },
  { "write_LP", impl_write_lp, TRUE },
  { "write_mps", impl_write_mps, TRUE },
  { "write_MPS", impl_write_mps, TRUE },
  { "write_params", impl_write_params, TRUE },
  { "write_XLI", impl_write_XLI, TRUE },

/* extra routines */

#if defined DEBUG || defined DEMO
  { "demo", impl_demo, FALSE },
#endif
  { "get_col_names", impl_get_col_name, TRUE },
  { "get_constr_types", impl_get_constr_type, TRUE },
  { "get_int", impl_is_int, TRUE },
  { "get_no_cols", impl_get_Ncolumns, TRUE },
  { "get_no_rows", impl_get_Nrows, TRUE },
  { "get_objective_name", impl_get_objective_name, TRUE },
  { "get_obj_fn", impl_get_obj_fn, TRUE },
  { "get_obj_fun", impl_get_obj_fn, TRUE },
  { "get_problem_name", impl_get_lp_name, TRUE },
  { "get_reduced_costs", impl_get_dual_solution, TRUE },
  { "get_row_names", impl_get_row_name, TRUE },
  { "mat_elm", impl_get_mat, TRUE },
  { "print_handle", impl_print_handle, FALSE },
  { "read_lp_file", impl_read_LP, FALSE },
  { "get_handle", impl_get_handle, FALSE },
  { "return_constants", impl_return_constants, FALSE },
};

#if defined WIN32
#  define signalconvention __cdecl
#else
#  define signalconvention
#endif

static void signalconvention SIGINT_func(int sig)
{
        interrupted = TRUE;
}

static void mainloop(structlpsolve *lpsolve)
{
	int	i;
        hashelem *hp;
        char verStr[128];

        interrupted = FALSE;
        signal(SIGINT, SIGINT_func);

        if (setjmp(lpsolve->lpsolvecaller.exit_mark) == 0) {
                char buf[bufsz];

                if (!initialized) {
        		/* Register the Exit Function */
                        registerExitFcn(lpsolve);

        		/* Allocate a string array to store command */

        		/* Allocate a string array to store error message */

                        /* create hashtable of all callbable commands to find them back quickly */
#if defined FORTIFY
                        Fortify_Disable(-1);
#endif
                        cmdhash = create_hash_table(sizeof(routines) / sizeof(*routines), 0);
                        for (i = 0; i < (int) (sizeof(routines)/sizeof(*routines)); i++)
                          	puthash(routines[i].cmd, i, NULL, cmdhash);

                        constanthash = create_hash_table(sizeof(constants) / sizeof(*constants), 0);
                        for (i = 0; i < (int) (sizeof(constants)/sizeof(*constants)); i++)
                          	puthash(constants[i].svalue, i, NULL, constanthash);

#if defined FORTIFY
                        Fortify_Disable(-2);
#endif

        		/* Initialise the lp array, and pointer to the last lp */

        		lp_last = -1;

                        if (!init_lpsolve_lib()) {
                                sprintf(buf, "Failed to initialise lpsolve library.%s", NEWLINE);
                                ErrMsgTxt(&lpsolve->lpsolvecaller, buf);
                        }

                        handlehash = NULL;

                	initialized = TRUE;
#                       if defined DEBUG
                        	Printf("Initialised%s", NEWLINE);
#                       endif
        	}

        	/* Get the first argument as a string matrix */

        	if (lpsolve->lpsolvecaller.nrhs < 1) {
        		int majorversion, minorversion, release, build;

        		/* ErrMsgTxt(&lpsolve->lpsolvecaller, "At least one command is required."); */
                        lp_solve_version(&majorversion, &minorversion, &release, &build);
                        if(lpsolve->lpsolvecaller.nlhs < 1) { //JONNY_EDIT
                            printSolverInfo();
                            //Printf(strdrivername "  " caller " Interface version " driverVERSION "%s" \
                            //       "using lpsolve version %d.%d.%d.%d%s%s" \
                            //       "Usage: ret = " strdrivername "(%sfunctionname%s, arg1, arg2, ...)%s",
                            //       NEWLINE, majorversion, minorversion, release, build, NEWLINE, NEWLINE, quotechar, quotechar, NEWLINE);
                        }
                        else {
                            sprintf(verStr,"%d.%d.%d.%d",majorversion,minorversion,release,build);
                            lpsolve->lpsolvecaller.plhs[0] = mxCreateString(verStr); 
                        }
                        return;
                }

        	GetString(&lpsolve->lpsolvecaller, NULL, 0, lpsolve->cmd, cmdsz, TRUE);

#               if defined DEBUG
               		Printf("%s%s", lpsolve->cmd, NEWLINE);
#               endif

        	/* Now call the required part of the lp toolkit */

                hp = findhash(lpsolve->cmd, cmdhash);
                if (hp == NULL) {
                        strcpy(buf, lpsolve->cmd);
                        strncat(buf,": Unimplemented.", bufsz);
                        ErrMsgTxt(&lpsolve->lpsolvecaller, buf);
                }
                i = hp->index;

                if (routines[i].needshandle) {
                        if (lpsolve->lpsolvecaller.nrhs < 2)
        		        ErrMsgTxt(&lpsolve->lpsolvecaller, "An lp handle is required.");

                        if (GetString(&lpsolve->lpsolvecaller, NULL, 1, buf, bufsz, FALSE)) {
                                if (handlehash != NULL)
                                	hp = findhash(buf, handlehash);
                                else
                                        hp = NULL;
                                if (hp == NULL) {
                                        char name[20 + bufsz];

                                        strcpy(name, buf);
                                        sprintf(buf, "Invalid model name: %s", name);
                                        ErrMsgTxt(&lpsolve->lpsolvecaller, buf);
                                }
                                lpsolve->h = hp->index;
                        }
                        else {
                                lpsolve->h = (int) GetRealScalar(&lpsolve->lpsolvecaller, 1);
                        }

                	if ((!handle_valid(lpsolve->h)) || ((lpsolve->lp = lp[lpsolve->h]) == NULL)) {
                        	strcpy(buf, lpsolve->cmd);
                        	strncat(buf, ": Invalid lp handle.", bufsz);
                		ErrMsgTxt(&lpsolve->lpsolvecaller, buf);
                	}
                }
                BEGIN_INTERRUPT_IMMEDIATELY_IN_FOREIGN_CODE;
                routines[i].routine(lpsolve);
                END_INTERRUPT_IMMEDIATELY_IN_FOREIGN_CODE;
        }
}

callerPrototype(drivername)
{
        structlpsolve lpsolve;
        struct structallocatedmemory *allocatedmemory;

        publicargs(&lpsolve);
        lpsolve.allocatedmemory = NULL;

        mainloop(&lpsolve);

        for (allocatedmemory = lpsolve.allocatedmemory; allocatedmemory != NULL;) {
                struct structallocatedmemory *allocatedmemory1;

                allocatedmemory1 = allocatedmemory;
                allocatedmemory = allocatedmemory->next;
                matFree(allocatedmemory1->ptr);
                matFree(allocatedmemory1);
        }

        ExitcallerPrototype(&lpsolve);
}

//Print Solver Information
void printSolverInfo()
{    
    int majorversion, minorversion, release, build;
    lp_solve_version(&majorversion, &minorversion, &release, &build);
    char vbuf[6]; getVSVer(vbuf); 
    
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" LP_SOLVE: Mixed Integer Linear Programming Solver [%d.%d.%d.%d, Built %s, VS%s]\n",majorversion,minorversion,release,build,__DATE__,vbuf);              
    mexPrintf("  - Released under the GNU Lesser General Public License: http://lpsolve.sourceforge.net/5.5/LGPL.htm\n");
    mexPrintf("  - Source available from: http://lpsolve.sourceforge.net/5.5/index.htm\n");
    
    mexPrintf("\n MEX Interface Peter Notebaert [Modified by J.Currie 2013]\n");
    mexPrintf("-----------------------------------------------------------\n");
}
