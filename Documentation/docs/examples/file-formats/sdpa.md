---
title: "SDPA and SeDuMi Models"
slug: "/examples/file-formats/sdpa/"
---

From OPTI v1.80 you can now read and write SDPA and SeDuMi models, in both sparse and dense format. SDPA and SeDuMi are the standard formats for storing semidefinite problems.

## SDPA Models
The SDPA file format was created as part of the   [SDPA Solver](http://sdpa.sourceforge.net/). It has definitions for both sparse and dense formats, however most models now use the sparse format. A good description of the SDPA format comes with the SDPA solver documentation, or the CSDP user's guide.

The  [SDPLIB](http://euler.nmt.edu/~brian/sdplib/sdplib.html) library is a collection of test semidefinite programs, some of which are included in OPTI (Test Problems/DAT-S). This collection uses the SDPA sparse format exclusively (as opposed to the dense format), and we will use a few of the problems below.

### Example 1: Loading a Sparse SDPA Model
To load a SDPA problem into OPTI, use the utility function shown below:

```matlab
>> prob = sdpRead('arch0.dat-s')
```

Returned will be a general MATLAB structure that you are free to process yourself, or it can be passed directly to opti:

```matlab
% Build an OPTI object of the returned problem 
OptSDPA = opti(prob)

% Solve the resulting OPTI object
[x,fval] = solve(OptSDPA)
```

You can load models from the MATLAB path, or if you supply an absolute path, to anywhere on your computer. The OPTI `sdpRead` function will return the problem in standard primal form, as is the OPTI convention.

### Example 2: Loading a Dense SDPA Model
A less used format is the SDPA dense format, where all elements (including zeros) are included. If you specify the file extension as just '.dat', OPTI will automatically read the input file in dense format. Alternatively you can specify the type manually, as below:

```matlab
>> probd = sdpRead('sdpademo','SDPA',1)
```

The above example shows the use of a couple of extra arguments. The second argument specifies the file type. This will be used if your file does not contain the extension within the first argument. The final argument indicates verbose mode, which will print information to the MATLAB command window as the file is processed. 

### Example 3: Writing an SDPA Model
Writing SDPA models is just as easy as reading them:

```matlab
>> sdpWrite(probd,'mysdpademo.dat-s')
```

Assuming you have been copying and pasting the examples so far, the above will write the first SDP problem back to a sparse SDPA file. You can also optionally write an OPTI object directly:

```matlab
>> write(OptSDPA,'myarch0.dat-s')
```

Which will write the OPTI object from Example 1 to a sparse SDPA file. Both methods will write to the current MATLAB directory, unless an absolute path is given.

## SeDuMi Models
The SeDuMi model is simply a model created in MATLAB and the model variables saved to a MATLAB .mat file. The  [SeDuMi Solver](http://sedumi.ie.lehigh.edu/) is a robust semidefinite programming solver which supports much more functionality than is currently interfaced by OPTI. 

### Example 1: Loading a SeDuMi Model
To load a SeDuMi problem into OPTI, use the utility function shown below:

```matlab
>> prob = sdpRead('sdp_arch0.mat')
```

Returned will be a general MATLAB structure that you are free to process yourself, or it can be passed directly to opti:

```matlab
% Build an OPTI object of the returned problem 
OptSDMI = opti(prob)

% Solve the resulting OPTI object
[x,fval] = solve(OptSDMI)
```

As before, you can load models from the MATLAB path, or if you supply an absolute path, to anywhere on your computer. <small>**However**</small> there is a major difference when loading SeDuMi models: the model will not be converted to OPTI standard primal form *unless* a solver other than SeDuMi is used. This way SeDuMi can be used to solve SeDuMi models without an intermediate conversion. You can also write SeDuMi problems using `sdpWrite` and specifying a .mat extension.

## Summary
While the methods shown here are attempt to read and write all SDPA and SeDuMi files, there may be errors in my code which mean the conversion is not correct. If you find a bug, feel free to report it on the Q&A Forum and I will try to fix it ASAP. Couple points to note, not all SeDuMi constraints are currently supported, and do not try to write large, sparse SDPA problems to dense files!
