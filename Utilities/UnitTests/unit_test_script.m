%% Individual Unit Test Running
% J.Currie October 2017
clc
clear

%% MKLJAC
clc
testMklJac = mklJac_tests;
res = run(testMklJac)

%% RMATHLIB
clc
testRmathlib = rmathlib_tests;
res = run(testRmathlib)