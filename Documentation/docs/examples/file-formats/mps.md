---
title: "MPS and LP Models"
slug: "/examples/file-formats/mps/"
---

Common formats for storing linear (and recently quadratic) models in text files are via the .mps and .lp formats. OPTI uses the COIN-OR CoinUtils package to be able to read and write both of these formats.

## MPS and QPS Models
Wikipedia gives a nice  [overview](http://en.wikipedia.org/wiki/MPS_%28format%29)  of the MPS specification, as does the  [lpsolve](http://lpsolve.sourceforge.net/5.5/mps-format.htm)  website. MPS is actually quite an old format, dating back to the days of punch cards! However it is still a useful format for linear and quadratic problems.

OPTI comes with a number of test mps files in the Test Problems/MPS/ directory you can use to look at it. It is a format that can be read by a human directly, although this is tedious!

### Example 1: Loading a MPS Model
To load an MPS problem into OPTI, use the utility function shown below:

```matlab
>> prob = coinRead('testLP.mps')
```

Returned will be a general MATLAB structure that you are free to process yourself, or it can be passed directly to opti:

```matlab
% Build an OPTI object of the returned problem 
OptMPS = opti(prob)

% Solve the resulting OPTI object
[x,fval] = solve(OptMPS)
```

You can load models from the MATLAB path, or if you supply an absolute path, to anywhere on your computer. The OPTI `coinRead` function will return the problem in row format (see the [constraint information](../../guides/advanced/cons.md) page), as is the default for COIN-OR solvers.

### Example 2: Loading a QPS Model
A QPS model is a quadratic program that is described in MPS format. In fact OPTI doesn't mind whether the file extension is .mps or .qps, if it contains a quadratic objective it will load it! An example is shown below:

```matlab
>> probq = coinRead('testQP2','qps',1)
```

The above example shows the use of a couple of extra arguments. The second argument specifies the file extension. This will be used if your file does not contain the extension within the first argument. The final argument indicates verbose mode, which will print information to the MATLAB command window as the file is processed. This is extremely useful for finding errors in MPS files!

Note in the above problem that the MPS default is that all lower bounds are 0, unless specified otherwise.

### Example 3: Loading a MILP with SOS
The MPS format also contains directives for integer variables, and Special Ordered Sets (SOS). The OPTI `coinRead` function will automatically read these and convert them to OPTI format, easy!

```matlab
>> Opt = opti(coinRead('testSOS2.mps'))
```

### Example 4: Writing a MPS Model
Writing MPS models is just as easy as reading them:

```matlab
>> coinWrite(prob,'mytest.mps')
```

Assuming you have been copying and pasting the examples so far, the above will write the first LP problem back to an MPS file. You can also optionally write an OPTI object directly:

```matlab
>> write(Opt,'mysostest.mps')
```

Which will write the OPTI object from Example 3 to an MPS file. Both methods will write to the current MATLAB directory, unless an absolute path is given.

## LP Models
The LP format accepted by OPTI is similar to that accepted by CPLEX, however only a very small subset of commands is recognised. Unlike the MPS format, an LP model is written as a set of equations which makes it much easier to read by eye.

Note the LP format appears to have multiple specifications, such as used by  [QSopt](http://www2.isye.gatech.edu/~wcook/qsopt/hlp/ff_lp_format.htm) and  [lpsolve](http://lpsolve.sourceforge.net/5.0/lp-format.htm). I believe the underlying parser only accepts models based on the  [CPLEX LP](http://lpsolve.sourceforge.net/5.0/CPLEX-format.htm) format.

### Example 5: Loading a LP Model
As in Example 1, to load an LP problem into OPTI, use the utility function shown below:

```matlab
>> prob = coinRead('prod.lp')
```

And returned will be a MATLAB structure containing the problem definition. The above function will determine what type of problem you are solving based on the file extension, or it can be passed as a second argument. As before, to solve the problem:

```matlab
% Build an OPTI object of the returned problem 
Opt = opti(prob)

% Solve the resulting OPTI object
[x,fval] = solve(Opt)
```

## Summary
While the methods shown here are attempt to read and write all MPS and LP files, there are a number of variations to the format which may not be accepted by the OPTI routines. However as the reading and writing is all done in compiled C++ and is super quick, no harm in trying!
