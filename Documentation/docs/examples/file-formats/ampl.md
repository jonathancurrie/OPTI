---
title: "AMPL .NL Models"
slug: "/examples/file-formats/ampl/"
---

AMPL has provided an interface called the AMPL Solver Library (ASL) which allows any model written in  [AMPL](http://www.ampl.com/) to be read and solved using any solver. 

Do note however OPTI cannot read an AMPL model directly. An intermediate file (.nl) must be generated from your AMPL model using the AMPL engine, which is then read by the ASL to allow a solver to solve it. This process is automatically run by AMPL every time you solve a problem, and from OPTI v1.77, is also automatically performed if the AMPL executable (ampl.exe) is on the MATLAB path.

Therefore you must still have a valid AMPL license (or use the free, problem size limited,  [student version](http://www.ampl.com/DOWNLOADS/details.html#WinStd)) to generate the .nl file from your model. OPTI then uses the ASL to read and subsequently solve your problem.

## Example 1: Generating an AMPL .NL file
Generating the required .nl file is pretty simple. From a command window call the following (assuming you are in the base directly of the ampl executable, or it is on the system path):

```dos
  ampl -og(1) (2) (3) (4)
```

Where the arguments are as follows:

1. Output file name
1. Input model .mod file
1. Optional data .dat file
1. Optional options .opt file

For example:
```dos
  ampl -ogdiet MODELS\diet.mod MODELS\diet.dat
```

Will generate the required .nl file (diet.nl) from the AMPL model diet.mod and dataset diet.dat. 

## Example 2: Loading an AMPL LP
To load an AMPL problem into OPTI, use the utility function shown below:

```matlab
prob = amplRead('diet.nl')
```

The function `amplRead` will automatically determine what type of problem is described in the model (from LP, MILP, QP, MIQP, QCQP, MIQCQP, NLP and MINLP) and create a return structure. You are free to process this structure yourself, or it can be passed directly to opti:

```matlab
% Build an OPTI object of the returned problem 
Opt = opti(prob)

% Solve the resulting OPTI object
[x,fval] = solve(Opt)
```

You can load models from the MATLAB path, or if you supply an absolute path, to anywhere on your computer. The OPTI `amplRead` function will return the problem in row format (see the [constraint information](../../guides/advanced/cons.md) page).

## Example 3: Loading an AMPL NLP
Unlike MPS or LP models, the AMPL NL format can also describe nonlinear models as well. And solving them is just as easy!

```matlab
% Load an AMPL NLP
prob = amplRead('ch3.nl')

% Build an OPTI object of the returned problem 
Opt = opti(prob)

% Solve the resulting OPTI object
[x,fval] = solve(Opt)

% Close ASL interface
asl('close')
```

A <small>**large**</small> benefit of using AMPL and the ASL is that all first and second derivatives are automatically calculated, and supplied to OPTI. This can greatly aid both the speed and robustness of the solver.

The final statement closes the ASL interace and is optional. When MATLAB closes all ASL interfaces you have open (via `amplRead`) will be automatically closed, releasing any used memory. However if you remember, it does not harm to include it in your code! This statement is only required when loading nonlinear models, as the objective and other callbacks remain part of the ASL.

## Example 4: Loading an AMPL MIP
AMPL can also accommodate Mixed Integer Programs (MIPs) and these can be read and solved by OPTI, both in linear and nonlinear forms:

```matlab
% Load an AMPL MILP
prob = amplRead('multmip1.nl')

% Build an OPTI object of the returned problem 
Opt = opti(prob)

% Solve the resulting OPTI object
[x,fval] = solve(Opt)
```

Note however AMPL will reorder the variables from the original model based on it's internal rules. Consult the AMPL PDFs supplied with OPTI for details on the reordering.

## Example 5: Loading an AMPL .mod file directly
As of OPTI v1.77 you can now directly solve AMPL .mod files from MATLAB, provided you have a licensed version of AMPL on the MATLAB path (i.e. ampl.exe). Thank you to Robert Fourer (ampl.com) for assisting with this development!

The following example is a non-convex MIQCQP which attempts to minimize losses in the paper industry. Using `amplRead` you can pass the model and data files as the first two arguments, then solve it as per the above examples:

```matlab
% Load an AMPL MIQCQP from Model File (but read as MINLP)
prob = amplRead('trimlon.mod','trimlon2.dat',[],1)

% Build an OPTI object of the returned problem 
Opt = opti(prob)

% Solve the resulting OPTI object
[x,fval] = solve(Opt)
```

Note the '1' passed as the fourth argument to `amplRead`. This new syntax from OPTI v1.79 tells the ASL interface to skip automatic identification of linear and quadratic problems, and instead treat the problem as purely nonlinear. This is required in this problem as OPTI may report an error if a quadratic constraint is not positive semidefinite (as it is in this case), and a global QCQP solver is not selected.

Remember the above example will only work if you have AMPL installed on your PC, and visible on the MATLAB path.

## Example 6: Using SCIP to solve an AMPL model {#scip-ampl}
If you wish to solve an AMPL model using SCIP, the existing OPTI AMPL interface would not provide an algebratic description, which is required by SCIP. Therefore from OPTI v1.78 I have added the native AMPL reader into the SCIP build. 

The new interface will be automatically selected by OPTI when solving an AMPL model, when SCIP is chosen as the solver, as shown below:

```matlab
% Load an AMPL MIQCQP from Model File
prob = amplRead('trimlon.mod','trimlon2.dat')

% Set SCIP as the solver
opts = optiset('solver','scip','display','iter');

% Build an OPTI object of the returned problem 
Opt = opti(prob,opts)

% Solve the resulting OPTI object
[x,fval] = solve(Opt)
```

The advantage here is that now we can solve non-convex nonlinear optimization problems and obtain a globally optimal solution. Remember SCIP does not solve problems with trigonometric functions. As with the above example, AMPL must be licensed and installed on your PC.

## Example 7: Extracting Constraint Linearity from an AMPL Model {#conlin}
From OPTI v1.79 `amplRead` now returns a new field, `conlin`, which contains a vector indicating the linearity of each constraint in your AMPL model. Using the Hock & Schittkowski #100 model as an example, the following code will identify linear, quadratic and nonlinear constraints:

```matlab
% Load an AMPL NLP from Model File
prob = amplRead('hs100.nl')

% Examine constraint linearity
prob.conlin
```

The above code will print the vector [-1,5,0]', where each element represents one of the three constraints. A nonlinear constraint will be indicated by < 0, a quadratic constraint by > 0 (technically the number of nz in the QC Q) and a linear constraint by 0. Remember AMPL may order the constraints differently than entered in the original .mod model.

## Summary
The AMPL interface is provided for people who have existing models in AMPL, but wish to try out OPTI or one of OPTI's solvers. While I have spent some time getting this interface going, there could still be some problems with it. If you find it is not returning the correct result (or crashing) I would be happy to take a look, if you can send me your model!
