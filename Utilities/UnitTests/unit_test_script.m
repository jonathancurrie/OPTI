%% Individual Unit Test Running
% J.Currie October 2017
clc
clear

%% OPTI Object
clc
testOPTI = opti_tests;
res = run(testOPTI)

%% MKLJAC
clc
testMklJac = mklJac_tests;
res = run(testMklJac)

%% RMATHLIB
clc
testRmathlib = rmathlib_tests;
res = run(testRmathlib)

%% CONFIDENCE LIMITS
clc
testConfidence = confidence_tests;
res = run(testConfidence)