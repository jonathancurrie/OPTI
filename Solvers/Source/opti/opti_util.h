#include <string.h>

void getVSVer(char vbuf[6])
{
    #if (_MSC_VER == 1500)
        strcpy(vbuf,"2008");
    #elif (_MSC_VER == 1600)
        strcpy(vbuf,"2010");
    #elif (_MSC_VER == 1700)
        strcpy(vbuf,"2012");
    #elif (_MSC_VER == 1800)
        strcpy(vbuf,"2013");
    #elif (_MSC_VER == 1900)
        strcpy(vbuf,"2015");
    #else
        strcpy(vbuf,"?");
    #endif
}