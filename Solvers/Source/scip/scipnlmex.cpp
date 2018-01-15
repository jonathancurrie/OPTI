/* SCIPMEX - A MATLAB MEX Interface to SCIP
 * Released Under the BSD 3-Clause License:
 * https://www.inverseproblem.co.nz/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2018
 * www.inverseproblem.co.nz
 */

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scipmex.h"
#include "mex.h"

//Function Prototypes
void debugPrintState(int state, int *args, int savexp, int savvar, int savpro, int varcnt);

//Message Array
extern char msgbuf[];

//Enumerations
enum {EMPTY=-2,READ,NUM,VAR,EXP,MUL,DIV,ADD,SUB,SQUARE,SQRT,POW,EXPNT,LOG,SIN,COS,TAN,MIN,MAX,ABS,SIGN};
const int OTHER = 98, EXIT = 99;

//Add Nonlinear Constraint to Problem
double addNonlinearCon(SCIP* scip, SCIP_VAR** vars, double *instr, size_t no_instr, double lhs, double rhs, double *x0, size_t nlno, bool isObj)
{
    //Internal Args
    SCIP_EXPR **expvars = NULL;     //expression array to represent variables
    SCIP_EXPR **exp = NULL;         //expression array, to incrementally build up full expression
    SCIP_EXPRTREE *exprtree = NULL; //resulting expression tree
    SCIP_CONS *nlcon;               //resulting constraint
    SCIP_VAR* nlobj;                //variable representing nonlinear objective (if used)
    int state = READ;               //state machine sate
    int index = 0;                  //variable index (picked up from instruction)
    int *varind;                    //array of all variable indicies
    int *unqind;                    //array of unique indicies
    int *expind;                    //array of expression variable indicies, matching to actual variables
    int varcnt = -1;                //keeps track of which expression variable we're up to
    int expno = 0;                  //index within expression array we are up to
    size_t no_var = 0;              //number of variables declared in instruction list
    size_t no_ops = 0;              //number of operations declared in instruction list
    size_t no_unq = 0;              //number of unique indicies found
    int args[2] = {EMPTY,EMPTY};    //used to keep track of previous arguments    
    int savexp = 0;                 //used to keep track of last expression before entering exp (op) expr
    int savexplst[MAX_DEPTH];       //list of saved expressions to process later
    int psavexplst = -1;            //position in list   
    int savvarlst[MAX_DEPTH];       //list of saved variables to process later (may not be needed...?)
    int psavvarlst = -1;            //position in list
    int vari = varcnt;              //varible index when evaluating functions   
    int savprolst[MAX_DEPTH*2];     //list of exp or var to process
    int psavprolst = -1;            //position in list
    size_t i, j;
    double num = 0;                 //to save constants passed
    double fval = 0;                //evaluation value
    double *lx0 = NULL;             //local x0 with just variables present in expression copied (and correct order)
    double one = 1.0, mone = -1.0, zero = 0.0;
    bool isunq = true;
    
    //Initialize Lists
    for(i = 0; i < MAX_DEPTH; i++) {
        savexplst[i] = EMPTY;
        savvarlst[i] = EMPTY;
    }
    for(i = 0; i < MAX_DEPTH*2; i++)
        savprolst[i] = EMPTY;
    
    //TEST CASE       
    //i.e. 2*x2*x4*x2, varind = [2 4 2], unqind = [2 4], var0 = x2, var1 = x4, expind = [0 1 0]
//     double instr[14] = {NUM,2.0,VAR,0,MUL,0,VAR,1,MUL,0,VAR,0,MUL,0};
//     double lhs = 1.0;
//     double rhs = 1.0;
//     no_instr = 14;

    //Determine number of variables and operations 
    //(note we have to have a new variable 'expression' every time a variable appears, even if its repeated)
    state = READ;
    for(i = 0; i < no_instr; i++) {
        switch(state)
        {
            case READ:
                if(instr[i] > EXP)
                    no_ops++;
                if(instr[i] == VAR)
                    no_var++;
                state = OTHER;
                break;
                
            case OTHER: //i.e. arg / misc
                state = READ;
                break;
        }
    }
    //Create Variable Index Vector
    varind = (int*)mxCalloc(no_var,sizeof(int));
    //Fill with indicies from instruction list
    j = 0; state = READ;
    for(i = 0; i < no_instr; i++) {
        switch(state)
        {
            case READ:
                if(instr[i] == VAR)
                    state = VAR;
                else
                    state = OTHER;
                break;
              
            case VAR:
                varind[j++] = (int)instr[i];
                state = READ;
                break;
                
            case OTHER:
                state = READ;
                break;
        }
    }
    //Find unique indicies to determine unique variables contained within the expression, for adding to the tree at the end
    unqind = (int*)mxCalloc(no_var,sizeof(int)); //make enough room for all to be unqiue
    //Fill in unique
    for(i = 0; i < no_var; i++) {
        if(!no_unq) { //no entries yet
            unqind[i] = varind[i];
            no_unq++;
        }
        else {
            //see if entry already exists
            isunq = true;
            for(j = 0; j < no_unq; j++) {
                if(unqind[j] == varind[i]) {
                    isunq = false;
                    break;
                }
            }
            //If is unique, add it
            if(isunq) 
                unqind[no_unq++] = varind[i];
        }        
    }
    //Allocate unique indicies a variable number (in this case just 0 to n), then allocate to expression index list
    expind = (int*)mxCalloc(no_var,sizeof(int));
    //Fill in unique indicies (to use in SCIPexprCreate() for variables)
    for(i = 0; i < no_var; i++) {
        for(j = 0; j < no_unq; j++) {
            if(varind[i] == unqind[j]) {
                expind[i] = j;
                break;
            }
        }
    }
    #ifdef DEBUG
        mexPrintf("\n---------------------------------------\nProcessing Nonlinear Expression\n---------------------------------------\n");
        mexPrintf("novar: %d, nounq: %d; no_ops: %d\n",no_var,no_unq,no_ops);

        //Print what we have found so far
        for(i = 0; i < no_var; i++) {
            mexPrintf("varind[%d] = %d\n",i,varind[i]);
        }
        mexPrintf("\n");
        for(i = 0; i < no_unq; i++) {
            mexPrintf("unqind[%d] = %d\n",i,unqind[i]);
        }
        mexPrintf("\n");
        for(i = 0; i < no_var; i++) {
            mexPrintf("expind[%d] = %d\n",i,expind[i]);
        }
        mexPrintf("\n");
    #endif
    
    //Create Expression Variable Memory
    SCIP_ERR( SCIPallocMemoryArray(scip,&expvars,no_var), "Error allocating expression variable memory");
    //Create Expression Memory
    SCIP_ERR( SCIPallocMemoryArray(scip,&exp,MAX(no_ops,1)), "Error allocating expression memory");
    //For each expression variable, create and index
    for(i = 0; i < no_var; i++)
        SCIPexprCreate(SCIPblkmem(scip), &expvars[i], SCIP_EXPR_VARIDX, expind[i]); //note we may have repeated indicies, as detailed above 
    
    //If objective, create an unbounded variable to add to objective, representing the nonlinear part
    if(isObj) 
    {
        SCIP_ERR( SCIPcreateVarBasic(scip, &nlobj, "nlobj", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS), "Error adding nonlinear objective variable");
        SCIP_ERR( SCIPaddVar(scip, nlobj), "Error adding nonlinear objective variable");
    }
    
    //Check for constant number objective
    if(no_instr == 2 && instr[0] == NUM) {
        SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[0], SCIP_EXPR_CONST, instr[1]), "Error creating constant objective / constraint expression"); 
        no_instr = 0; //skip below
        #ifdef DEBUG
            mexPrintf("Found constant objective / constraint, skipping expression tree building\n");
        #endif
    }
    //Check for single variable objective / constraint
    if(no_instr == 2 && instr[0] == VAR) {
        SCIP_ERR( SCIPexprCreateLinear(SCIPblkmem(scip), &exp[0], 1, &expvars[0], &one, 0.0), "Error creating linear expression of a single variable");
        no_instr = 0; //skip below
        #ifdef DEBUG
            mexPrintf("Found single variable objective / constraint, skipping expression tree building\n");
        #endif
    }
    
    //Begin processing instruction list
    state = READ;
    for( i = 0; i < no_instr; i++)
    {        
        #ifdef DEBUG
            debugPrintState(state, args, psavexplst, psavvarlst, psavprolst, varcnt);
        #endif
        switch(state)
        {
            case READ:  //Read Instruction
                state = (int)instr[i];
                break;
                
            case NUM:   //Number passed, collect value
                num     = instr[i];
                state   = READ;
                if(args[0] == EMPTY) //all empty
                    args[0] = NUM;
                else if(args[1] == EMPTY) //tag onto end
                    args[1] = NUM;
                else if(args[0] == EXP) { //all full, save expression for use later
                    if(psavexplst >= MAX_DEPTH)
                        mexErrMsgTxt("Maximum function depth exceeded [expression list]");
                    //Save expression number in list
                    savexplst[++psavexplst] = expno-1;
                    //Update process list
                    savprolst[++psavprolst] = EXP;
                    //Correct args
                    args[0] = args[1];
                    args[1] = NUM;                    
                }
                else if(args[0] == VAR) { //all full, save variable for use later
                    if(psavvarlst >= MAX_DEPTH)
                        mexErrMsgTxt("Maximum function depth exceeded [variable list]");
                    //Save variable number in list
                    savvarlst[++psavvarlst] = varcnt-1;
                    //Update process list
                    savprolst[++psavprolst] = VAR;
                    //Correct args
                    args[0] = args[1];
                    args[1] = NUM;
                }
                else
                    mexErrMsgTxt("Error in order of instructions");
                break;     
                
            case VAR:   //Variable index passed, increment our variable counter
                varcnt++;
                state   = READ;
                if(args[0] == EMPTY) //all empty
                    args[0] = VAR;
                else if(args[1] == EMPTY) //tag onto end
                    args[1] = VAR;
                else if(args[0] == EXP) { //all full
                    if(psavexplst >= MAX_DEPTH)
                        mexErrMsgTxt("Maximum function depth exceeded [expression list]");
                    //Save expression number in list
                    savexplst[++psavexplst] = expno-1;
                    //Update process list
                    savprolst[++psavprolst] = EXP;
                    args[0] = args[1];
                    args[1] = VAR;
                }
                else if(args[0] == VAR) { //all full but with variables e.g x2 - x1^2
                    if(psavvarlst >= MAX_DEPTH)
                        mexErrMsgTxt("Maximum function depth exceeded [variable list]");
                    //Save variable number in list
                    if(instr[i-3] == VAR && instr[i-5] == VAR) //check last two instructions left hand col (this will be at least the third instr)
                        savvarlst[++psavvarlst] = varcnt-2;
                    else
                        savvarlst[++psavvarlst] = varcnt-1;
                    //Update process list
                    savprolst[++psavprolst] = VAR;
                    //Correct args
                    args[0] = args[1];
                    args[1] = VAR;
                }
                else
                    mexErrMsgTxt("Error in order of instructions");
                break;
            
            //All two operand operators
            case MUL:   
            case DIV:
            case ADD:
            case SUB:
            case POW:
                
                //Read variable index
                vari = varcnt;                
               
                //Check for waiting variable or expression to process
                if(args[1] == EMPTY && psavprolst >= 0) {
                    //Read if the last record is a var or exp, then check if we actually have one to process
                    int toProcess = savprolst[psavprolst--];
                    if(toProcess == VAR) {
                        if(psavvarlst < 0) //none to process, error
                            mexErrMsgTxt("Error reading variable to process - Process List indicates variable to process, but not present in variable list");
                        //Collect variable number to process
                        vari = savvarlst[psavvarlst];
                        if(vari == EMPTY || vari < 0)
                            mexErrMsgTxt("Error processing waiting variable, found empty or negative index");  
                        #ifdef DEBUG
                            mexPrintf("-----\nProcessing Waiting Variable: %d [Index %d, %d remaining]\n-----\n",vari,varind[vari],psavvarlst);
                        #endif
                        //Correct args
                        args[1] = VAR;
                        instr[i] = 1; //will flip if subtraction or division  
                        psavvarlst--;
                    }
                    else if(toProcess == EXP) {
                        if(psavexplst < 0) //none to process, error
                            mexErrMsgTxt("Error reading expression to process - Process List indicates expression to process, but not present in expression list");
                        #ifdef DEBUG
                            mexPrintf("-----\nProcessing Waiting Expression: %d [%d remaining]\n-----\n",savexplst[psavexplst],psavexplst);
                        #endif
                        args[1] = EXP; 
                    }
                    else
                        mexErrMsgTxt("Unknown entry in process list");
                }
                
                //Check we have two arguments
                if(args[0] == EMPTY || args[1] == EMPTY)
                    mexErrMsgTxt("Error attempting to create nonlinear expression, operator doesn't have two operands");

                //NUM (op) VAR (linear except division)
                if(args[0] == NUM && args[1] == VAR)
                { 
                    switch(state)
                    {
                        case MUL: SCIP_ERR( SCIPexprCreateLinear(SCIPblkmem(scip), &exp[expno], 1, &expvars[vari], &num, 0.0), "Error creating linear multiply expression (num * var)"); break;
                        case ADD: SCIP_ERR( SCIPexprCreateLinear(SCIPblkmem(scip), &exp[expno], 1, &expvars[vari], &one, num), "Error creating linear add expression (num + var)"); break;
                        case SUB: SCIP_ERR( SCIPexprCreateLinear(SCIPblkmem(scip), &exp[expno], 1, &expvars[vari], &mone, num), "Error creating linear subtract expression (num - var)"); break;
                        case DIV: 
                            SCIP_EXPR *divexp; //intermediate expression
                            SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &divexp, SCIP_EXPR_CONST, num), "Error creating division intermediate expression"); 
                            SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_DIV, divexp, expvars[vari]), "Error creating divide expression (num / var)"); break;
                        case POW:
                            mexErrMsgTxt("You cannot use POWER with the exponent as a variable. For x^y use exp(y*log(x))");
                        default:
                            mexErrMsgTxt("Operator not implemented yet for NUM (op) VAR!");
                    }
                    args[0] = EXP;
                    args[1] = EMPTY;
                }
                //VAR (op) NUM (linear except pow)
                else if(args[0] == VAR && args[1] == NUM)
                { 
                    switch(state)
                    {
                        case MUL: SCIP_ERR( SCIPexprCreateLinear(SCIPblkmem(scip), &exp[expno], 1, &expvars[vari], &num, 0.0), "Error creating linear multiply expression (var * num)"); break;
                        case ADD: SCIP_ERR( SCIPexprCreateLinear(SCIPblkmem(scip), &exp[expno], 1, &expvars[vari], &one, num), "Error creating linear add expression (var + num)"); break;
                        case SUB: SCIP_ERR( SCIPexprCreateLinear(SCIPblkmem(scip), &exp[expno], 1, &expvars[vari], &one, -num), "Error creating linear subtract expression (var - num)"); break;
                        case DIV: SCIP_ERR( SCIPexprCreateLinear(SCIPblkmem(scip), &exp[expno], 1, &expvars[vari], &(num = 1/num), 0.0), "Error creating linear divide expression (var * 1/num)"); break;
                        case POW:
                            //Check for integer power (may need to think of a better way)
                            if((double)((int)num) == num) {
                                if(num==2) { //squared
                                    SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_SQUARE, expvars[vari], (int)num), "Error creating squared expression");}
                                else {
                                    SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_INTPOWER, expvars[vari], (int)num), "Error creating integer power expression (var ^ inum)");}
                            }
                            else {
                                /*if(num < 0) //seems to work?
                                    mexErrMsgTxt("This interface does not support a variable raised to a negative real power");
                                else*/ 
                                    SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_REALPOWER, expvars[vari], num), "Error creating power expression (var ^ num)");
                            }                               
                            break;
                            
                        default:
                            mexErrMsgTxt("Operator not implemented yet for VAR (op) NUM!");
                    }
                    args[0] = EXP;
                    args[1] = EMPTY;
                }
                //VAR (op) VAR
                else if(args[0] == VAR && args[1] == VAR)
                {
                    switch(state)
                    {
                        case MUL: SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_MUL, expvars[vari-1], expvars[vari]), "Error creating mul expression (var * var)"); break;
                        case DIV: SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_DIV, expvars[vari-1], expvars[vari]), "Error creating div expression (var / var)"); break;
                        case ADD: SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_PLUS, expvars[vari-1], expvars[vari]), "Error creating add expression (var + var)"); break;
                        case SUB: SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_MINUS, expvars[vari-1], expvars[vari]), "Error creating sub expression (var - var)"); break;
                        case POW:
                            mexErrMsgTxt("You cannot use POWER with the exponent as a variable. For x^y use exp(y*log(x))");    
                        default:
                            mexErrMsgTxt("Operator not implemented yet for VAR (op) VAR!");
                    }
                    args[0] = EXP;
                    args[1] = EMPTY;
                }
                //EXP (op) NUM
                else if(args[0] == EXP && args[1] == NUM)
                { 
                    switch(state)
                    {
                        case MUL: SCIP_ERR( SCIPexprCreateLinear(SCIPblkmem(scip), &exp[expno], 1, &exp[expno-1], &num, 0.0), "Error creating linear multiply expression (exp * num)"); break;                        
                        case ADD: SCIP_ERR( SCIPexprCreateLinear(SCIPblkmem(scip), &exp[expno], 1, &exp[expno-1], &one, num), "Error creating linear add expression (exp + num)"); break;
                        case SUB: 
                            //Check for flipped args
                            if(instr[i] == 1) { //flip
                                SCIP_EXPR *divexp; //intermediate expression
                                SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &divexp, SCIP_EXPR_CONST, num), "Error creating subtraction intermediate expression"); 
                                SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_MINUS, divexp, exp[expno-1]), "Error creating subtract expression (num - exp)");
                            }
                            else
                                SCIP_ERR( SCIPexprCreateLinear(SCIPblkmem(scip), &exp[expno], 1, &exp[expno-1], &one, -num), "Error creating linear subtract expression (exp - num)");
                            break;                            
                        case DIV: 
                            //Check for flipped args
                            if(instr[i] == 1) { //flip
                                SCIP_EXPR *divexp; //intermediate expression
                                SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &divexp, SCIP_EXPR_CONST, num), "Error creating division intermediate expression"); 
                                SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_DIV, divexp, exp[expno-1]), "Error creating divide expression (num / exp)");
                            }
                            else
                                SCIP_ERR( SCIPexprCreateLinear(SCIPblkmem(scip), &exp[expno], 1, &exp[expno-1], &(num = 1/num), 0.0), "Error creating linear divide expression (exp * 1/num)"); 
                            break;
                        case POW:
                            //Check for integer power (may need to think of a better way)
                            if((double)((int)num) == num) {
                                SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_INTPOWER, exp[expno-1], (int)num), "Error creating integer power expression (exp ^ inum)");}
                            else {
                                /*if(num < 0) //seems to work?
                                    mexErrMsgTxt("This interface does not support an expression raised to a negative real power");
                                else*/ 
                                    SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_REALPOWER, exp[expno-1], num), "Error creating power expression (exp ^ num)");
                            }                               
                            break;
                            
                        default:
                            mexErrMsgTxt("Operator not implemented yet for NUM (op) VAR!");
                    }
                    args[0] = EXP;
                    args[1] = EMPTY;
                }
                //EXP (op) VAR
                else if(args[0] == EXP && args[1] == VAR) 
                {                  
                    switch(state)
                    {
                        case MUL: SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_MUL, exp[expno-1], expvars[vari]), "Error creating mul expression (exp * var)"); break;                        
                        case ADD: SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_PLUS, exp[expno-1], expvars[vari]), "Error creating add expression (exp + var)"); break;
                        case SUB: 
                            //Check for flipped args
                            if(instr[i] == 1) {  //flip
                                SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_MINUS, expvars[vari], exp[expno-1]), "Error creating subtract expression (var - exp)");}
                            else
                                SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_MINUS, exp[expno-1], expvars[vari]), "Error creating subtract expression (exp - var)");
                            break; 
                        case DIV:
                            //Check for flipped args
                            if(instr[i] == 1) { //flip
                                SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_DIV, expvars[vari], exp[expno-1]), "Error creating div expression (var / exp)");}
                            else
                                SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_DIV, exp[expno-1], expvars[vari]), "Error creating div expression (exp / var)"); 
                            break;
                        case POW:
                            mexErrMsgTxt("You cannot use POWER with the exponent as a variable. For x^y use exp(y*log(x))");
                            
                        default:
                            mexErrMsgTxt("Operator not implemented yet for EXP (op) VAR!");
                    }
                    args[0] = EXP;
                    args[1] = EMPTY;
                }
                //EXP (op) EXP
                else if(args[0] == EXP && args[1] == EXP)
                {
                    //Collect expression number to process
                    savexp = savexplst[psavexplst];
                    if(savexp == EMPTY || savexp < 0)
                        mexErrMsgTxt("Error processing waiting expression, found empty index");
                    
                    switch(state)
                    {
                        case ADD: SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_PLUS, exp[savexp], exp[expno-1]), "Error creating add expression (exp + exp)"); break;
                        case SUB: 
                            //Check for flipped args
                            if(instr[i] == 1) { //flip
                                SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_MINUS, exp[expno-1], exp[savexp]), "Error creating sub expression (exp_old - exp_new)");}
                            else
                                SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_MINUS, exp[savexp], exp[expno-1]), "Error creating sub expression (exp_new - exp_old)"); 
                            break;
                        case DIV:
                            //Check for flipped args
                            if(instr[i] == 1) { //flip
                                SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_DIV, exp[expno-1], exp[savexp]), "Error creating div expression (exp_old / exp_new)");}
                            else
                                SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_DIV, exp[savexp], exp[expno-1]), "Error creating div expression (exp_new / exp_old)"); 
                            break;
                        case MUL: SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_MUL, exp[savexp], exp[expno-1]), "Error creating add expression (exp * exp)"); break;
                            
                        default:
                            mexErrMsgTxt("Unexpected operator for combining expressions, currently only +, -, * and / supported");
                    }                    
                    args[0] = EXP;
                    args[1] = EMPTY;
                    psavexplst--;
                }
                else
                    mexErrMsgTxt("Grouping of arguments not implemented yet");
                
                if(i == (no_instr-1))
                    state = EXIT;
                else {
                    state = READ;
                    expno++; //increment for next operator
                }                
                break;
                
            //All single operand functions
            case SQUARE:
            case SQRT:
            case EXPNT:
            case LOG:
            case ABS:
                
                //Read variable index
                vari = varcnt; 
                
                //Check we have one argument
                if(args[0] == EMPTY)
                    mexErrMsgTxt("Error attempting to create nonlinear expression, function doesn't have an operand");
                //Check for unprocessed exp [e.g. (3-x(2)).*log(x(1))]
                if(args[0] == EXP && args[1] == VAR)
                {
                    //Save expression to process later
                    if(psavexplst >= MAX_DEPTH)
                        mexErrMsgTxt("Maximum function depth exceeded [expression list]");
                    //Save expression number in list
                    savexplst[++psavexplst] = expno-1;
                    //Update process list
                    savprolst[++psavprolst] = EXP;
                    //Correct args
                    args[0] = VAR;
                    args[1] = EMPTY;
                }
                //Check for unprocessed var [e.g. x2 - exp(x1)]
                if(args[0] == VAR && args[1] == VAR)
                {   
                    //Save variable to process later
                    if(psavvarlst >= MAX_DEPTH)
                        mexErrMsgTxt("Maximum function depth exceeded [variable list]");
                    //Save variable number in list
                    savvarlst[++psavvarlst] = varcnt-1;
                    //Update process list
                    savprolst[++psavprolst] = VAR;
                    //Continue processing function, we will process this variable next free round
                }                

                //FCN ( VAR ) enum {EMPTY=-2,READ,NUM,VAR,EXP,MUL,DIV,ADD,SUB,SQUARE,SQRT,POW,EXPNT,LOG,SIN,COS,TAN,MIN,MAX,ABS,SIGN};
                if(args[0] == VAR)
                { 
                    //Check variable index realistic
                    if(vari < 0)
                        mexErrMsgTxt("Variable index is negative when creating fcn(var)");

                    switch(state)
                    {
                        case SQRT:  SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_SQRT, expvars[vari]), "Error creating sqrt expression sqrt(var)"); break;
                        case EXPNT: SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_EXP, expvars[vari]), "Error creating exponential expression exp(var)"); break;
                        case LOG:   SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_LOG, expvars[vari]), "Error creating natural log expression log(var)"); break;
                        case ABS:   SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_ABS, expvars[vari]), "Error creating absolute value expression abs(var)"); break;                       
                        default:
                            mexErrMsgTxt("Operator not implemented yet for FCN ( VAR )!");
                    }
                    args[0] = EXP;
                    args[1] = EMPTY;
                }
                //FCN ( EXP )
                else if(args[0] == EXP)
                {
                    switch(state)
                    {
                        case SQRT:  SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_SQRT, exp[expno-1]), "Error creating sqrt expression sqrt(exp)"); break;
                        case EXPNT: SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_EXP, exp[expno-1]), "Error creating exponential expression exp(exp)"); break;
                        case LOG:   SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_LOG, exp[expno-1]), "Error creating natural log expression log(exp)"); break;
                        case ABS:   SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp[expno], SCIP_EXPR_ABS, exp[expno-1]), "Error creating absolute value expression abs(exp)"); break;                       
                        default:
                            mexErrMsgTxt("Operator not implemented yet for FCN ( EXP )!");
                    }
                    args[0] = EXP;
                    args[1] = EMPTY;
                }
                else
                    mexErrMsgTxt("Unknown argument to operate function on, only options are VAR (variable) or EXP (expression)");
                
                if(i == (no_instr-1))
                    state = EXIT;
                else {
                    state = READ;
                    expno++; //increment for next operator
                }
                break;
                
            case MIN: 
            case MAX:
                mexErrMsgTxt("Max and Min not currently implemented (in this interface)");
            case SIN:
            case COS:
            case TAN:
                mexErrMsgTxt("Trigonometric functions not currently implemented (in SCIP)");
            case SIGN:
                mexErrMsgTxt("Sign not currently implemented (in SCIP)");
                
            case EXIT:
                break; //exit switch case
                
            default:
                mexErrMsgTxt("Unknown (or out of order) instruction");
        }       
    }
    #ifdef DEBUG
        debugPrintState(state, args, psavexplst, psavvarlst, psavprolst, varcnt);
        mexPrintf("\n---------------------------------------\nSummary at Expression Tree Create:\n");
        mexPrintf("expno:  %3d [should equal %3d]\nvarcnt: %3d [should equal %3d]\n",expno,no_ops-1,varcnt,no_var-1);
        mexPrintf("---------------------------------------\n");
    #endif
    //Create Expression Tree to hold our final expression
    SCIP_ERR( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, exp[expno], no_unq, 0, NULL), "Error creating expression tree");
    //Create Variable Array of variables to add (can't seem to find a method to does this one by one)
    SCIP_VAR **unqvars = NULL;
    SCIP_ERR( SCIPallocMemoryArray(scip,&unqvars,no_unq), "Error allocating unique variable array");
    for(i = 0; i < no_unq; i++)
        unqvars[i] = vars[unqind[i]]; //store address of variable into double pointer array
    //Set the variables within the tree
    SCIP_ERR( SCIPexprtreeSetVars(exprtree, no_unq, unqvars), "Error setting expression tree var");
    
    //'Validate' the constraint if initial guess supplied
    if(x0 != NULL) {
        //Copy in x0 variables used in the expression
        lx0 = (double*)mxCalloc(no_var,sizeof(double));
        //Note variables appear in the same order as found in expression
        for(i = 0; i < no_unq; i++)
            lx0[i] = x0[unqind[i]];
        SCIPexprtreeEval(exprtree,lx0,&fval);
        //Free memory
        mxFree(lx0);
    }
    
    //Create the Nonlinear Constraint, Add It, then Release It
    if(isObj)
        sprintf(msgbuf,"NonlinearObj%d",nlno);
    else
        sprintf(msgbuf,"NonlinearExp%d",nlno);
    SCIP_ERR( SCIPcreateConsBasicNonlinear(scip, &nlcon, msgbuf, 0, NULL, NULL, 1, &exprtree, &one, lhs, rhs), "Error creating nonlinear constraint!");  
    SCIPexprtreeFree(&exprtree);
    //If nonlinear objective, add to objective
    if(isObj)
        SCIP_ERR( SCIPaddLinearVarNonlinear(scip, nlcon, nlobj, -1.0), "Error adding nonlinear objective linear term");
    //Continue to free the constraint  
    SCIP_ERR( SCIPaddCons(scip, nlcon), "Error adding nonlinear constraint");
    SCIP_ERR( SCIPreleaseCons(scip, &nlcon), "Error freeing nonlinear constraint"); 
              
    //Memory clean up
    SCIPfreeMemoryArray(scip, &exp);
    SCIPfreeMemoryArray(scip, &expvars);   
    SCIPfreeMemoryArray(scip, &unqvars);
    mxFree(expind);
    mxFree(varind);
    mxFree(unqind);
    
    return fval;
}

void debugPrintState(int state, int *args, int savexp, int savvar, int savpro, int varcnt)
{
    char *strs[18] = {"EXIT","EMPTY","READ","NUM","VAR","EXPRSN","MUL","DIV","ADD","SUB","SQUARE","SQRT","POWER","EXPNT","LOG","SIN","COS","TAN"};
    if(state == 99) state = -3;
    mexPrintf("State: %-8s ARG 0: %-8s ARG 1: %-8s PSAVEXP: %3d PSAVVAR: %3d PSAVPRO: %3d VARCNT: %3d\n",strs[state+3],strs[args[0]+3],strs[args[1]+3],savexp,savvar,savpro,varcnt);
}

    //Test code
    /*SCIP_EXPR *exp0, *exp1, *var0, *var1;
    SCIPexprCreate(SCIPblkmem(scip), &var0, SCIP_EXPR_VARIDX, 0);
    SCIPexprCreate(SCIPblkmem(scip), &var1, SCIP_EXPR_VARIDX, 1);
    
    SCIP_ERR( SCIPexprCreateLinear(SCIPblkmem(scip), &exp0, 1, &var0, num, 0.0), "Error creating linear expression"); //2*x0        
    SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp1, SCIP_EXPR_MUL, exp0, var1), "Error creating mul expression"); //(2*x0)*x1
    
    SCIP_ERR( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, exp1, 2, 0, NULL), "Error creating expression tree");*/
  
        /*SCIP_EXPR *exp0, *exp1, *exp2, *exp3, *var0, *var1, *var2;
    SCIPexprCreate(SCIPblkmem(scip), &var0, SCIP_EXPR_VARIDX, 0);
    SCIPexprCreate(SCIPblkmem(scip), &var1, SCIP_EXPR_VARIDX, 1);
    SCIPexprCreate(SCIPblkmem(scip), &var2, SCIP_EXPR_VARIDX, 1);
    
    num = 2.5;
    SCIP_ERR( SCIPexprCreateLinear(SCIPblkmem(scip), &exp0, 1, &var0, &num, 0.0), "Error creating linear expression"); //2.5*x0        
    SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp1, SCIP_EXPR_MUL, exp0, var1), "Error creating mul expression"); //(2.5*x0)*x1
    
    num = 0.2;
    SCIP_ERR( SCIPexprCreateLinear(SCIPblkmem(scip), &exp2, 1, &var2, &num, 0.0), "Error creating linear expression"); //0.2*x1        
    SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp3, SCIP_EXPR_MINUS, exp1, exp2), "Error creating mul expression"); //((2.5*x0)*x1) - (0.2*x1)
    
    SCIP_ERR( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, exp3, 2, 0, NULL), "Error creating expression tree");
    SCIP_ERR( SCIPexprtreeSetVars(exprtree, 2, vars), "Error setting expression tree var");
    
        SCIP_EXPR *exp0, *exp1, *var0, *var1;
    num = 0.2;
    SCIPexprCreate(SCIPblkmem(scip), &var0, SCIP_EXPR_VARIDX, 0);
    SCIPexprCreate(SCIPblkmem(scip), &var1, SCIP_EXPR_CONST, num);
    
    num = 2.0;
    SCIP_ERR( SCIPexprCreateLinear(SCIPblkmem(scip), &exp0, 1, &var0, &num, 0.0), "Error creating linear expression"); //2*x0        
    SCIP_ERR( SCIPexprCreate(SCIPblkmem(scip), &exp1, SCIP_EXPR_DIV, var1, exp0), "Error creating div expression"); //0.2 / (2*x0)
    
    SCIP_ERR( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, exp1, 2, 0, NULL), "Error creating expression tree");
    SCIP_ERR( SCIPexprtreeSetVars(exprtree, 2, vars), "Error setting expression tree var");*/
