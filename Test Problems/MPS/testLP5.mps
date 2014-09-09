NAME          TESTPROB
ROWS
 N  COST
 L  1
 G  2
 E  3
COLUMNS
    XONE      COST                 1        1                 1
    XONE      2                    1
    YTWO      COST                 4        1                 1
    YTWO      3                    -1
    ZTHREE    COST                 9        2                 1
    ZTHREE    3                    1
RHS
    RHS1      1                    5        2                 10
    RHS1      3                    7
BOUNDS
 UP BND1      XONE                 4
 LO BND1      YTWO                -1
 UP BND1      YTWO                 1
ENDATA
