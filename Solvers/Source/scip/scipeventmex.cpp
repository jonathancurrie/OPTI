/* SCIPMEX - A MATLAB MEX Interface to SCIP
 * Released Under the BSD 3-Clause License:
 * https://www.inverseproblem.co.nz/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2018
 * www.inverseproblem.co.nz
 */

#include "scipmex.h"
#include "mex.h"
#include <signal.h>

//Ctrl-C Detection
#ifdef __cplusplus
    extern "C" bool utIsInterruptPending();
    extern "C" void utSetInterruptPending(bool);
#else
    extern bool utIsInterruptPending();
    extern void utSetInterruptPending(bool);
#endif

//Executed when adding the event
static SCIP_DECL_EVENTINIT(eventInitCtrlC)
{
    //Notify SCIP to add event
    SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODESOLVED, eventhdlr, NULL, NULL) );
    //Return
    return SCIP_OKAY;   
}

//Executed when removing the event
static SCIP_DECL_EVENTEXIT(eventExitCtrlC)
{
    //Notify SCIP to drop event
    SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_NODESOLVED, eventhdlr, NULL, -1) );
    //Return
    return SCIP_OKAY;
}

//Executed when event occurs
static SCIP_DECL_EVENTEXEC(eventExecCtrlC)
{
    //Check for Ctrl-C
    if (utIsInterruptPending()) {
        utSetInterruptPending(false); /* clear Ctrl-C status */
        mexPrintf("\nCtrl-C Detected. Exiting SCIP...\n\n");
        raise(SIGINT);
    }
    return SCIP_OKAY; //always OK - otherwise we don't get intermediate answer
}

//Called from scipmex to include our event
SCIP_RETCODE SCIPincludeCtrlCEventHdlr(SCIP* scip)
{
    SCIP_EVENTHDLR* eventhdlr = NULL;
    
    //Create Event Handler
    SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, "CtrlCMatlab", "Catching Ctrl-C From Matlab", eventExecCtrlC, NULL) );
    
    //Setup callbacks
    SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitCtrlC) );
    SCIP_CALL( SCIPsetEventhdlrExit(scip, eventhdlr, eventExitCtrlC) );
    
    return SCIP_OKAY;
}