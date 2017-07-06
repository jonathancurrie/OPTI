
#if (_MSC_VER == 1500)
    #define VS_VER "2008"
#elif (_MSC_VER == 1600)
    #define VS_VER "2010"
#elif (_MSC_VER == 1700)
    #define VS_VER "2012"
#elif (_MSC_VER == 1800)
    #define VS_VER "2013"
#elif (_MSC_VER == 1900)
    #define VS_VER "2015"
#elif (_MSC_VER == 1910)
    #define VS_VER "2017"
#else
    #define VS_VER "?"
#endif

#ifndef ML_VER
#error Please define the MATLAB version when building the MEX file
#endif
#ifndef OPTI_VER
#error Please define the OPTI version when building the MEX file
#endif

#define STR_EXPAND_OPTI(tok) #tok
#define stringify(tok) STR_EXPAND_OPTI(tok)
#define PRINT_BUILD_INFO mexPrintf("  - Built %s, Visual Studio %s, MATLAB %s, OPTI v%s\n", __DATE__, VS_VER, stringify(ML_VER), stringify(OPTI_VER));