---
title: "GAMS .GMS Models"
slug: "/examples/file-formats/gams/"
---

[GAMS](http://www.gams.com/) is a high-level algebraic modelling language and system that allows optimization using a range of advanced solvers. Using OPTI you can convert your MATLAB optimization problem into a GAMS model, allowing it to be solved under their modelling environment. This feature could be useful for debugging and validating solutions, or for trying different solvers via  [NEOS](http://www.neos-server.org/neos/solvers/index.html).

OPTI requires one of two possible packages for writing GAMS models. If you have SCIP in your OPTI distribution (i.e. you have downloaded the academic version) then you already have everything you need! SCIP does all the hard work, and you can turn any SCIP compatible linear, quadratic or nonlinear model into a GAMS model.

On the other hand if you do not have SCIP, you will need to download the latest  [GAMS distribution](http://www.gams.com/download/). This is free, and allows solving of small problems. It also includes the GAMS Data Exchange (GDX) routines for MATLAB, which OPTI will utilize for writing model data to. Simply download and install GAMS, then add the GAMS directory to the MATLAB path (e.g. C:\GAMS\win64\24.0). OPTI includes routines for writing problems up to MIQCQP complexity to GAMS models without SCIP using GDX. Nonlinear models cannot be written to a GAMS model without SCIP currently.

## Example 1: Writing a GAMS .GMS file
Regardless of whether you have SCIP, or GDX, or both, the function call is the same. The following example writes a simple linear program to a GAMS model:

```matlab
% Linear Objective & Constraints
f = -[6 5]';
A = ([1,4; 6,4; 2, -5]); 
b = [16;28;6];    

%Build OPTI Object
Opt = opti('grad',f,'ineq',A,b,'bounds',[0;0],[10;10]);

% Write to GAMS Model
write(Opt,'OPTI_LP1.gms')
```

If you have SCIP, the model will be contained with single .gms file. If you are using GDX, then as well as the .gms model, you will also have a opti_lp1.gdx file, which contains the model parameters.

## Example 2: Solving a Non-Convex QCQP via NEOS
Let's assume you don't have SCIP in your OPTI distribution, yet you want to try it out on a non-convex QCQP. First, let's enter the model:

```matlab
% Quadratic Objective
H = eye(2);                
f = -[2 2]';       
         
% Linear Constraints
A = [-1,1; 1,3];            
b = [2;5];    
lb = [0;0];     
           
% Quadratic Constraints (qrl <= x'Qx + l'x <= qru)
Q = {[1 0; 0 1]            
      [1 0; 0 1]};
l = {[0;-2]; [-2;2]};
qrl = {3; 1};          %QC1 is double sided, QC2 is an equality
qru = {5; 1};
```

Now if you try and build this model without specifying a solver, OPTI will tell you it is non-convex and throw an error. This is because you don't have a native QCQP solver that can solve this problem. However if you specify IPOPT as the solver, OPTI will do a little bit of work to transform the problem, and allow you to create the object:

```matlab
% Create OPTI Object
Opt = opti('qp',H,f,'ineq',A,b,'lb',lb,'qcrow',Q,l,qrl,qru,...
           'options',optiset('solver','ipopt'))

% Write the QCQP Problem
write(Opt,'ncQCQP.gms')
```

You will now have two files in the current directory, ncQCQP.gms and ncQCQP.gdx. However to allow NEOS to accept the .gdx file and select SCIP as the solver, we will need to make two changes to ncQCQP.gms. 

- On line 23 (or thereabouts) will be "`$gdxin ... \ncQCQP.gdx`". Change this line to "`$gdxin in.gdx`".
- Before the "`Solve`" statement, add the line "`OPTION QCP = SCIP;`"

Now we will solve our model using SCIP, hosted on NEOS [here](http://www.neos-server.org/neos/solvers/minco:scip/GAMS.html). Click on the previous link, select your model and GDX file, ensure you agree with the Terms and Conditions of NEOS, then click "Submit to NEOS". In a few seconds you will see the result is -2.817, which for OPTI SCIP users, will be the same.

If you do have SCIP installed, and are following this example, you will not have the GDX file so you will not need to make the first change (in.gdx), but you will need to modify the solve statement to your desired solver.

## Example 3: Writing an MINLP
For all general nonlinear problems (i.e. not quadratic) you will need SCIP to convert your OPTI problem into a GAMS model. This is because part of the SCIP interface is used to parse the MATLAB model, as well SCIP itself performs the model writing. Also note GAMS currently applies default bounds of 0-100 on integer variables, so watch out!

The following example shows it is just as easy (with SCIP of course) to write nonlinear problems to GAMS models:
```matlab
% Objective
obj = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;

% Linear Constraints
A = [-1 1]; b = -1;
Aeq = [1 1]; beq = 5; 
lb = [0;0]; ub = [4;4];

%Integer Constraints
xtype = 'IC';

% Create OPTI Object
Opt = opti('obj',obj,'bounds',lb,ub,'xtype',xtype,'ineq',A,b,'eq',Aeq,beq)

% Write the MINLP Problem
write(Opt,'minlp.gms')
```

## Reading GAMS Models
The next obvious question is - can I read GAMS models into MATLAB? The answer is sort-of. GAMS can be called directly from MATLAB (using the [gams](http://www.gams.com/dd/docs/solvers/convert.pdf) command), however I have not tried it. You can also convert a GAMS model to an AMPL model using the [GAMS Convert Utility](http://www.gams.com/dd/docs/solvers/convert.pdf), then read it into MATLAB using the [AMPL reading commands](./ampl.md) supplied with OPTI.

Basically it is difficult, and don't think many people would have a use for this feature. If you really want to be able to read GAMS models into MATLAB I can look at automating the above, simply send me an email.
