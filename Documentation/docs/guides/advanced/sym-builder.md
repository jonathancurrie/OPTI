---
title: "SymBuilder"
slug: "/guides/advanced/sym-builder/"
---

The idea of SymBuilder (Symbolic Builder) is to take an algebraic equation / symbolic orientated approach to building optimization problems in MATLAB. This page presents a brief introduction to SymBuilder and its core functionality.

The main benefits of using SymBuilder are:
- It will automatically identify what type of problem you are solving (LP/QP/NLP), regardless of how you enter the problem (i.e. no need for vectors/matrices for LP/QP, equations are fine).
- If solving a nonlinear problem, it will automatically generate full analytical first and second derivatives for use by the optimizer, including sparsity patterns.
- If full analytical derivatives are too expensive to evaluate, SymBuilder can generate a C++ model and automatically compile it against [CppAD](../../reference/utilities/cpp-ad.md) to generate automatic derivatives.
- Nonlinear models can be automatically converted to MATLAB code, C code or C++ code, including MEX file interfaces, where applicable.
- By maintaining a Symbolic version of your model, the problem is automatically symbolically simplified by the Symbolic Toolbox.

There are some drawbacks: it can be slow on large problems, it requires you to re-write your model and it requires the Symbolic Toolbox, however generally I find it performs well. 

<small>**Note**</small> If you get an error message similar to `"Undefined function 'symb_cb'..."` when using SymBuilder, then try type `rehash` at the MATLAB command line. This will force MATLAB to find any recently generated callback functions.

## Example 1: Creating a SymBuilder Object
The SymBuilder object is created with no input arguments:

```matlab
% Create SymBuilder Object
B = SymBuilder();
```

or if you want to control the verbosity level, with a single input boolean

```matlab
% Create SymBuilder Object with suppressed command line output 
Bs = SymBuilder(false);
```

where `false` suppresses all SymBuilder information messages, and `true` (the default) prints all messages.

Both objects created above are empty and ready to accept a problem definition. The object is designed so that the problem is entered with minimal pre-processing done as the problem is entered (for efficiency), and then finally when it is 'Built' (see the end of this example), it constructs the complete optimization problem.

### Adding an Objective
To add an objective to the model, use the `AddObj` method. The objective equation is supplied as a string, or a Symbolic Toolbox symbolic variable expression. 

```matlab
% Add as a String (scalar expression only)
B.AddObj('-6*x1 -5*x2')

% Add as a Symbolic Expression
x = sym('x',[2 1]); %define symbolic variables
Bs.AddObj(-[6;5]'*x)
```

By default SymBuilder always minimizes the objective function. If you supply the equation as a string you do not need to identify the variables, the method will automatically determine the optimization variables when the object is built.

Also note that the linear objective is just an example, the object can accept just about any objective function that the Symbolic Toolbox can understand.

### Adding Constraints
Constraints are added using the `AddCon` method. As with the objective, general constraints can be added as strings or symbolic expressions. When supplied as symbolic expressions, the constraint numerical upper and lower bounds (cl <= c(x) <= cu) must be supplied as the 2nd and 3rd arguments. 

```matlab
% Add as strings
B.AddCon('x1 + 4*x2 <= 16');
B.AddCon('6*x1 + 4*x2 <= 28')

% Add as Symbolic Expressions with specified cl, cu
Bs.AddCon(x(1) + 4*x(2),-Inf,16);
Bs.AddCon(6*x(1) + 4*x(2),-Inf,28)
```

Note that using both methods, double sided constraints can also be declared. For equality constraints, use `"="` for a string, and set `cl` = `cu` for symbolic expressions. It is also possible to add a vector symbolic expression.

As with the objective, the linear constraints above are just an example. Quadratic and nonlinear constraints are fine as well.

### Adding Rectangular Bounds
Rectangular bounds can be added in three ways: 1) individually per variable, 2) as a vectorized expression, or 3) as symbolic variables with individual lb and ub specified.

```matlab
% Method 1) Individually as Strings
B.AddBound('0 <= x1 <= 10');
B.AddBound('0 <= x2 <= 10')

% Method 2) Vectorized, 'x' is recognised as belonging to all 'xn' variables
B.AddBound('0 <= x <= 10');

% Method 3) As a symbolic variable vector with numerical bounds
Bs.AddBounds(x,[0;0],[10;10])
```

If you declare bounds on a variable multiple times, the last declaration will be used.

### Building the SymBuilder Object
There are two methods available for turning the templated object into an optimization problem; `Draft` which computes just 1st derivatives, and `Build` which computes both 1st and 2nd derivatives. `Build` is useful for small problems, or when you want to attempt to identify a quadratic problem. `Draft` is useful when the problem is linear (as in this one), or when generating Symbolic 2nd Derivatives is too expensive.

```matlab
% Build Optimization Problem
Build(B)
```

At this point SymBuilder has automatically determined the type of problem entered, and is ready to construct an (internal) OPTI object to solve the problem. The automatic identification step will identify all OPTI problems excluding NLS/SNLE/SDP.

### Solving the Symbuilder Object
Once built (or drafted), call the `Solve` method to automatically generate the internal OPTI object and solve the problem.

```matlab
% Solve SymBuilder Optimization Problem
[x,fval,ef,info] = Solve(B)
```

Note if solving a nonlinear problem you must supply `x0` as the second argument to `Solve`. For the above example, we are solving a simple linear program, thus `x0` is not required.

## Example 2: Mixed Integer Quadratic Program
As stated, SymBuilder does not just work for linear programs, it was in fact designed for NLPs, including integer variants. The following is a [MIQP](../../examples/problem-types/miqp.md), automatically constructed and solved with SymBuilder:

```matlab
% New SymBuilder Object
B = SymBuilder();

% Add Quadratic Objective
B.AddObj('0.5*x1^2 + 0.5*x2^2 + 0.5*x3^2 - 2*x1 - 3*x2 - x3');

% Add Linear Constraints
B.AddCon('x1 + x2 + x3 <= 1');
B.AddCon('3*x1 - 2*x2 - 3*x3 <= 1');
B.AddCon('x1 - 3*x2 + 2*x2 <= 1');

% Add Integer Constraint (also possible with Symbolic Variables)
B.AddInteger('x2 = I');

% Build Object
Build(B)

% Solve Problem
[x,fval,ef,info] = Solve(B)
```

Note SymBuilder correctly identified the problem as a MIQP, which we can see by examining the Built object:

```matlab
>> B

------------------------------------------------------
SymBuilder Object
 BUILT in 0.072s with:
 -  3 variables
 -  1 objective
      -  0 linear
      -  1 quadratic
      -  0 nonlinear
 -  3 constraint(s)
      -  3 linear
      -  0 quadratic
      -  0 nonlinear
 -  0 bound(s)
 -  1 integer variable(s)
      -  1 integer
      -  0 binary
------------------------------------------------------
```

If we were to call `Draft` instead of `Build`, SymBuilder cannot establish that this is a quadratic program (by not being able to examine the second derivatives), and instead will treat it as a [MINLP](../../examples/problem-types/minlp.md), which could be much less efficient. 

## Example 3: Adding Constants / Parameters
There may be instances when variables are actually constants (or parameters), but you don't want to hard-code the number into an equation string. To illustrate, consider the following MINLP:

```matlab
% New SymBuilder Object
B = SymBuilder();

% Add Nonlinear Objective with 2 Constants
B.AddObj('sin(pi*x1/a1)*cos(pi*x2/a2)');

% Add Linear Constraints
B.AddCon('-x1 + 2.5*x2 <= 1');
B.AddCon('x1 + 2.5*x2 <= -15');

% Add Integer Constraint (note also vectorized)
B.AddInteger('x = I');
```

In the above example, the objective function contains two constants, `a1` and `a2`. We know they are constants because we wrote the program, it is nothing to do with the choice of variable name! To tell SymBuilder that they are constants, and not variables (which it will assume by default), use the following commands:

```matlab
% Declare Constants with Numerical Values
B.AddConstant('a1',12);
B.AddConstant('a2',16);
```

If we now build the object, and examine the generated equations, we will see SymBuilder has substituted the constants into the objective function:

```matlab
% Build It
Build(B)

% Inspect resulting equations
B.sobj

>>  cos((pi*x2)/16)*sin((pi*x1)/12)
                        2.5*x2 - x1
                        x1 + 2.5*x2
```

This technique allows you to add problem dependent parameters into common model equations.

## Example 4: Adding Expressions
SymBuilder also allows you to use intermediate expressions when constructing objective or constraints. Using the same example problem as above, the below example combines both intermediate expressions and constants.

```matlab
% New SymBuilder Object
B = SymBuilder();

% Add Nonlinear Objective with 2 Intermediate Expressions
B.AddObj('sin(e1)*cos(e2)');

% Add Linear Constraints
B.AddCon('-x1 + 2.5*x2 <= 1');
B.AddCon('x1 + 2.5*x2 <= -15');

% Add Integer Constraint
B.AddInteger('x = I');

% Add Intemediate Variable Expressions (with Constants for this example)
B.AddExpression('e1 = pi*x1/a1');
B.AddExpression('e2 = pi*x2/a2');

% Define Constants
B.AddConstant('a1',12);
B.AddConstant('a2',16);

% Add Now Build It
Build(B)
```

This technique allows you to simplify long equations into a series of shorter, hopefully simpler equations.

## Example 5: Displaying Results {#dispres}
Inspecting long solution variable vectors can be tedious and error prone, so SymBuilder allows the user to create 'Result Groups', and then add variables or expressions to display to the user once the object is solved.

To illustrate, consider the following (very) hypothetical vehicle routing problem:

```matlab
% New SymBuilder Object (verbose=false)
B = SymBuilder(false);

% Add Quadratic Objective
B.AddObj('(x1-x2)^2 + (x2+x3-2)^2 + (x4-1)^2 + (x5-1)^2');

% Add Linear Constraints
B.AddCon('x1 + 3*x2 = 5');
B.AddCon('x3 + x4 - 2*x5 = 0');
B.AddCon('x2 - x5 = 0');
```

Obviously all our variables are called x1, x2, etc, however they represent more meaningful real quantities for this problem. To create variable groups and assign them, use `AddResultGroup` and `AddResultExp`:

```matlab
% Add Result Groups (the letter is arbitrary, any group name can be used)
B.AddResultGroup('A','Fuel Usage [l]');
B.AddResultGroup('B','Distance [km]');

% Add Result Expressions (Group:Name, Variable or Expression)
B.AddResultExp('A:Truck 1','x1');
B.AddResultExp('A:Truck 2','x2');
B.AddResultExp('B:Route 1','x3');
B.AddResultExp('B:Route 2','x4');
B.AddResultExp('B:Route 3','x5-x3'); %expressions with variables OK too
```

Now we will build the object, solve it, and call the `Results` method to automatically display formatted results:

```matlab
% Build Object
Build(B);

% Solve
Solve(B,[ 2.5 0.5 2 -1 0.5 ]);

% Display Formatted Results
Results(B)

------------------------------------------------------
SymBuilder Optimization Results
 SOLVED in 0.00017777s
 STATUS: Optimal
 ITERATIONS: 1

 COST: 0.255814
------------------------------------------------------
A: Fuel Usage [l]
  - Truck 1      = 1.4419       
  - Truck 2      = 1.186        

B: Distance [km]
  - Route 1      = 1.093        
  - Route 2      = 1.2791       
  - Route 3      = 0.093023     

------------------------------------------------------
```

## Example 6: Code Generation {#code-gen}
A new feature in OPTI v2.10 is the ability to generate C or C++ code from a *nonlinear* SymBuilder model. While this functionality is basically just using the Symbolic Toolbox `ccode` method, SymBuilder wraps it and automates the entire MEX file generation process.

To illustrate, consider the QCQP (treated for this example as an NLP) below:

```matlab
% Enter Problem
B = SymBuilder();
B.AddObj('(x1-x2)^2 + (x2+x3-2)^2 + (x4-1)^2 + (x5-1)^2');
B.AddCon('x1^2 + 3*x2 = 4');
B.AddCon('x3 + x4 - 2*x5 = 0');
B.AddCon('x2 - x5 = 0');

% Draft Model
Draft(B)

% Initial Guess
x0 = [ 2.5 0.5 2 -1 0.5 ];
```

### Generating C-Code with Analytical Derivatives
To generate a C-code MEX file with analytical derivatives from the above model, simply set the callback mode (`'cbmode'`) in `symbset` as follows:

```matlab
% Generate C-Code Model and Solve Problem
[x,fval,ef,info] = Solve(B,x0,symbset('cbmode','ccode'))
```

The above code will generate a C source file of your optimization problem, automatically invoke the MEX compiler to compile the source to a dynamic library, then create and solve an optimization problem with the compiled library. 

A few points to remember on the above:
- For large problems generating the C-code can take some time, but the execution of the code thereafter should be much faster than equivalent MATLAB code (of course always problem dependent).
- You must use Visual Studio or Windows SDK for compiling the file, LCC does not work.

### Generating C++ Code with Automatic Derivatives
A neat new feature in SymBuilder is to exploit AD of a C++ version of your model. OPTI v2.10 includes [CppAD](../../reference/utilities/cpp-ad.md) as part of the distribution, and will automatically compile it into your model to generate all derivatives (sparse by default), including patterns. Simply change the callback mode:

```matlab
% Generate C++ Code Model with CppAD and Solve Problem
[x,fval,ef,info] = Solve(B,x0,symbset('cbmode','cppad'))
```

and SymBuilder will generate a C++ version of your model, automatically compile it against CppAD, then solve the optimization problem. For problems with particularly complex (dense) derivatives, AD can be an attractive option. Also note there is no point calling `Build` on a model you are going to solve with CppAD, as it internally generates all derivatives (unless of course your model is actually quadratic, in which case CppAD would not be used). 

A few points to remember on the above:
- Compiling large problems against CppAD can take a *long* time. This is due to CppAD being a template library, and the C++ optimizer appears to have a hard time optimizing the generated code. I am aware of this and working on a fix with the author of CppAD.
- You must use Visual Studio or Intel C++ to compile CppAD models, LCC and Windows SDK do not work.

### What Does SymBuilder do by Default?
By default, SymBuilder's `cbmode` is on `auto`. This tries to make an intelligent guess at what would be the most effective callback mode, between `mcode` and `cppad`. For small or particularly sparse problems, `mcode` is chosen, while for larger dense problems, `cppad` is used (unless a suitable compiler is not available). 

Oh - and if you're looking for the generated C/C++ source, SymBuilder automatically deletes it. If you would like to keep it, change the `symbset` option `'srckeep'` to `'yes'`.

## Example 7: Extracting Problem Data
If you want to extract the functions or data generated by SymBuilder, simply use `GetLinProb` for LPs, `GetQuadProb` for QP/QCQPs or `GetNLProb` for NLPs to return an `optiprob` structure with the problem data.

```matlab
% New SymBuilder Object
B = SymBuilder(false);

% Add Nonlinear Objective
B.AddObj('sin(pi*x1/12)*cos(pi*x2/16)');

% Add Linear Constraints
B.AddCon('-x1 + 2.5*x2 <= 1');
B.AddCon('x1 + 2.5*x2 <= -15');

% Build Object
Build(B);

% Get Nonlinear Problem Data
nlprob = GetNLProb(B)
```

## Example 8: SymBuilder Options
SymBuilder uses the `symbset` routine to control problem generation and solving. Use this routine to generate an options structure suitable for supplying to `Solve` or `GetNLProb`.

```matlab
% Print all Options
symbset

% Set Solver, Display to 'iter'
sopts = symbset('solver','clp','display','iter');
```

## Example 9: `optisym` {#optisym}
`optisym` is a simple API to allow users familiar with `fmincon` type functions to use SymBuilder and leverage its functionality. Simply pass normal MATLAB functions to `optisym` and it internally converts them to Symbolic Toolbox expressions and then to a SymBuilder object.

`optisym` has the calling form:
```matlab
[optiOpj,SymBobj] = optisym(fun,x0,lb,ub,con,cl,cu,xtype,sopts,verbose)
```

As with normal SymBuilder models, if you pass a LP (even as function handle), `optisym` will identify the problem as such and solve it as a LP. This provides a convenient method to simplifying suitable models and solving them as efficiently as possible.

Note as per standard SymBuilder functions, this API is only compatible with functions that the Symbolic Toolbox can understand and parse.

To illustrate, consider NLP Hock & Schittkowski #71

```matlab
% Objective Function (min fun(x))
fun = @(x) x(1)*x(4)*sum(x(1:3)) + x(3);

% Nonlinear Constraint Function (cl <= con(x) <= cu)
con = @(x) [ prod(x);
             sum(x.^2)];
cl = [25;40];
cu = [Inf;40];

% Bounds + x0
lb = ones(4,1);
ub = 5*ones(4,1);
x0 = [1 5 5 1]';

% Build OPTI Object
Opt = optisym(fun,x0,lb,ub,con,cl,cu)

% Solve
[x,fval,exitflag,info] = solve(Opt)
```

For more examples, see `test_probs_sym.m` in Test Problems/Development.

## Example 10: `Draft` vs `Build`
In most of the previous examples we used the `Build` command to generate our SymBuilder model. `Build` takes the symbolic problem and generates both symbolic first and second derivatives, which can be very time consuming on large or dense problems. The generation of second derivatives can help an optimizer, but it may just not be worth the time waiting to generate them, and thus an approximation may be quicker (e.g. IPOPT can use an internal L-BFGS update).

To skip generating second derivatives, you can use the `Draft` command instead. To illustrate, consider the following QCQP: 

```matlab
B = SymBuilder();

% Enter Objective
B.AddObj('(x1-x2)^2 + (x2+x3-2)^2 + (x4-1)^2 + (x5-1)^2');

% Add Constraints
B.AddCon('x1^2 + 3*x2 = 4');
B.AddCon('x3 + x4 - 2*x5 = 0');
B.AddCon('x2 - x5 = 0');
```

Now if we use `Draft` we get:

```matlab
Draft(B)

------------------------------------------------------
SymBuilder Object
 DRAFT in 0.023s with:
 -  5 variables
 -  1 objective
      -  0 linear
      -  1 nonlinear
 -  3 constraint(s)
      -  2 linear
      -  1 nonlinear
 -  0 bound(s)
 -  0 integer variables(s)
------------------------------------------------------
```

versus `Build`:
```matlab
Build(B)

------------------------------------------------------
SymBuilder Object
 BUILT in 0.095s with:
 -  5 variables
 -  1 objective
      -  0 linear
      -  1 quadratic
      -  0 nonlinear
 -  3 constraint(s)
      -  2 linear
      -  1 quadratic
      -  0 nonlinear
 -  0 bound(s)
 -  0 integer variables(s)
------------------------------------------------------
```

Immediately we can see due to the generation of second derivatives, SymBuilder was able to determine we are solving a quadratically constrained quadratic program, rather than a general nonlinear program as found by `Draft`. This can be a substantial performance increase when solved as a QCQP! Therefore if you suspect your problem is actually quadratic, `Build` can render a more efficient problem definition, and one that does not require the generation of a callback function (the problem is stored as a collection of matrices and vectors).

However when using OPTI, NLP solvers are normally used for solving QCQPs, thus as described above, the time waiting for the extra second derivative information may not be worthwhile. Let's assume you just called `Draft`, and therefore want to solve it as a NLP:

```matlab
% Enter Problem
B = SymBuilder();
B.AddObj('(x1-x2)^2 + (x2+x3-2)^2 + (x4-1)^2 + (x5-1)^2');
B.AddCon('x1^2 + 3*x2 = 4');
B.AddCon('x3 + x4 - 2*x5 = 0');
B.AddCon('x2 - x5 = 0');
% Draft
Draft(B)

% Solve
x0 = [ 2.5 0.5 2 -1 0.5 ];
[x,fval,ef,info] = Solve(B,x0)
```

You will receive a warning indicating 2nd derivatives were not found and therefore not added to the callback function. This is because the default callback function mode is to include 2nd derivatives. To disable this warning, disable the need for 2nd derivatives:

```matlab
% Solve without generating 2nd derivative callback
[x,fval,ef,info] = Solve(B,x0,symbset('use2ndDerivs','no'))
```

This example shows that you can control how the derivatives are generated and used with a SymBuilder model. Normally I use `Build`, however as shown in the next example, `Draft` can be useful if using a global solver (which doesn't need derivatives) or a callback with automatic differentiation.

## Summary
While the SymBuilder syntax is limited, being able to automatically generate exact first and second derivatives is quite attractive, especially on real problems. If generating the callback functions for nonlinear problems is taking too long, you can experiment with the various code-generation modes to see if an improvement can be found.

As an example how I used it in my work in steam utility optimization, I created a derived class from the SymBuilder object:

```matlab
classdef SymUtility < SymBuilder
%SYMUTILITY  Create a SYMUTILITY object for creating symbolic utility system optimization
```

then added a series of custom methods which expanded the SymBuilder object for modelling common steam system equipment, e.g.

```matlab
function [con1,con2] = AddWHB(B,JStm,names,duty,bfwH,stmT,stmP,bdr,eff)
%ADDWHB Add Waste Heat Boiler to Symbolic Builder Model [Duty Based]
%
%   AddWHB(JStm,names,duty,bfwH,stmT,stmP,bdr,eff)
%
%   names are {Unit Name,msteam,mbfw,h}
%   eff is optional, defaults to 0.9

    if(nargin < 9), eff = 0.9; end
    %Thermo Calcs
    WHB_H = JStm.HPT(stmP,stmT);
    WHB_BDH = JStm.HPX(stmP,0);
    WHB_eff = eff; %assume fixed efficiency
    %Energy Balance
    con1 = sprintf('%1.15g - %s = 0',(duty*WHB_eff)/...
                   (WHB_H + bdr*WHB_BDH - bfwH*(bdr + 1))*3.6,names{2});
    %Mass Balance
    con2 = sprintf('%s-(1+%1.15g)*%s = 0',names{3},bdr,names{2}); 
    %Add To Model
    B.AddConstant(names{4},WHB_H);
    B.AddConstant([names{1} '_Q'],duty);
    B.AddConstant([names{1} '_Eff'],eff);
    B.AddCon(con1); B.AddCon(con2);
    B.AddResultExp(['F:' names{1}],names{2}); 
end
```

As shown above, you can use `sprintf` to generate equations which can then be added to the SymBuilder object. The result of this approach is that I could create my optimization problem as:

```matlab
% Waste Heat Boiler
U.AddWHB(JStm,{'WHB','m13','m37','WHB_H'},WHB_Q,BFW_H,WHB_T,MP_P,0.03,WHB_Eff);
% Turbo Generators
U.AddTurboGen(JStm,{'TG1','m4','h3','btg1'},HP_P,LP_P,TG1_QMax,HP_H);
% Back Pressure Turbines
U.AddBPT(JStm,{'BT2','m6','h3','bt2'},HP_P,LP_P,BT2_Q,BT2_Eff,HP_H);
```

then finally at the end, call `Build`, then `Solve`, then `Results`, and easily inspect my optimal solution!
