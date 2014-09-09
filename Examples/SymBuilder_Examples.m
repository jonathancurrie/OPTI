%% Symbolic Builder (SymBuilder) Examples
%
%   SymBuilder is used for developing algebraic / symbolic optimization
%   problems, with an emphasis on automatically developing analytical
%   derivative expressions. It is also a basic modelling language,
%   extending the functionality of the Symbolic Toolbox (required).
%
%   Copyright (C) 2014 Jonathan Currie (I2C2)
clear

% There is also a page on the Wiki which supplements this example:
web('http://i2c2.aut.ac.nz/Wiki/OPTI/index.php/Advanced/SymBuilder');

% NOTE - If you get an error similar to "Undefined function 'symb_cb'" when
% using SymBuilder then try type "rehash" at the MATLAB command line. This 
% will force MATLAB to find any recently generated callback functions.

%% Example 1: Creating a SymBuilder Object
% The SymBuilder object is created with no arguments (default
% configuration), or a single boolean which specifies the verbosity level
% (true = print messages = default).

clc
B = SymBuilder();       %Prints Building Messages
Bs = SymBuilder(false); %Supresses Text Output

% Both objects created above are empty, and ready to accept a problem.
display(B)

%% Adding An Objective
% To add an objective to the model, use the AddObj() method. The objective
% equation is supplied as a string, or a symbolic variable expression. Note
% when supplying as strings, variables do not need to be identified prior
% to use in the object - SymBuilder does this automatically.

clc
% Add as a String (scalar expression only)
B.AddObj('-6*x1 -5*x2')

% Add as a Symbolic Expression
x = sym('x',[2 1]);
Bs.AddObj(-[6;5]'*x)

%% Adding Constraints
% Constraints are added using the AddCon() method. As with the objective,
% general constraints can be added as strings or symbolic expressions. When
% supplied as symbolic expressions, the constraint upper and lower bounds
% (cl <= c(x) <= cu) must be supplied as the 2nd and 3rd arguments. Note
% double sided constraints can also be used.

clc
% Add as strings
B.AddCon('x1 + 4*x2 <= 16');
B.AddCon('6*x1 + 4*x2 <= 28')

% Add as Symbolic Expressions with specified cl, cu
Bs.AddCon(x(1) + 4*x(2),-Inf,16);
Bs.AddCon(6*x(1) + 4*x(2),-Inf,28)

%% Adding Rectangular Bounds
% Rectangular bounds can be added in three ways: 1) individually per
% variable, 2) as a vectorized expression, or 3) as symbolic variables with
% individual lb and ub specified.

clc
% Method 1) Individually as Strings
B.AddBound('0 <= x1 <= 10');
B.AddBound('0 <= x2 <= 10')

% Method 2) Vectorized, 'x' is recognised as belonging to all xn
B.AddBound('0 <= x <= 10');

% Method 3) As a symbolic variable vector with numerical bounds
Bs.AddBounds(x,[0;0],[10;10])

%% Building the SymBuilder Object
% There are two methods available for turning the templated object into an
% optimization problem; Draft() which computes just 1st derivatives, and
% Build() which computes both 1st and 2nd derivatives. Build() is useful
% for small problems, or when you want to attempt to identify a quadratic
% problem. Draft() is useful when the problem is linear (as in this one),
% or when generating Symbolic 2nd Derivatives is too expensive.

clc
% Object thus far
B
% Build Optimization Problem
Build(B)

% At this point SymBuilder has automatically determined the type of problem
% entered, and is ready to construct an OPTI object to solve the problem.
% The automatic identification step will identify all OPTI problems
% excluding NLS/SNLE/SDP.

%% Solving the SymBuilder Object
% Once built (or drafted), call the Solve() method to automatically
% generate a OPTI object and solve the problem.

clc
[x,fval,ef,info] = Solve(B)


%% Example 2: Mixed Integer Quadratic Program Example
% SymBuilder does not just work for linear programs, it was in fact designed
% for QPs and NLPs, including integer variants. The following example is
% a MIQP, automatically constructed and solved with SymBuilder:

clc
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

% Note SymBuilder correctly identified the problem as MIQP and used an
% appropriate solver to solve it.


%% Example 3: Adding Constants
% There may be instances when variables are actually constants, but you
% don't want to hard-code the number into an equation string. To
% illustrate, consider the following example MINLP:

clc
% New SymBuilder Object
B = SymBuilder();
% Add Nonlinear Objective with 2 Constants
B.AddObj('sin(pi*x1/a1)*cos(pi*x2/a2)');
% Add Linear Constraints
B.AddCon('-x1 + 2.5*x2 <= 1');
B.AddCon('x1 + 2.5*x2 <= -15');
% Add Integer Constraint
B.AddInteger('x = I');

% If we build the object now, we will see that problem entered has 4
% variables. However let us define the values for constants/parameters a1
% and a2:
B.AddConstant('a1',12);
B.AddConstant('a2',16);

% Add Now Build It
Build(B)

% Inspect resulting equations
B.sobj

% This technique allows you to add problem dependent parameters into common
% equations.

%% Example 4: Adding Expressions
% SymBuilder also allows you to use intermediate expressions when
% constructing objective or constraints. Using the same example problem as 
% above, the below example combines both intermediate expressions and
% constants.

clc
% New SymBuilder Object
B = SymBuilder();
% Add Nonlinear Objective with 2 Intermediate Expressions
B.AddObj('sin(e1)*cos(e2)');
% Add Linear Constraints
B.AddCon('-x1 + 2.5*x2 <= 1');
B.AddCon('x1 + 2.5*x2 <= -15');
% Add Integer Constraint
B.AddInteger('x = I');

% Add Intemediate Variable Expressions with Constants
B.AddExpression('e1 = pi*x1/a1');
B.AddExpression('e2 = pi*x2/a2');
% Define Constants
B.AddConstant('a1',12);
B.AddConstant('a2',16);

% Add Now Build It
Build(B)

% Inspect resulting equations
B.sobj

% This technique allows you to simplify long equations into a series of
% shorter, simpler equations.

%% Example 5: Displaying Results
% Inspecting long solution variable vectors can be tedious and error prone,
% so SymBuilder allows the user to create Result Groups, and then add
% variables or expressions to display to the user. 

clc
% New SymBuilder Object (verbose=false)
B = SymBuilder(false);
% Add Quadratic Objective
B.AddObj('(x1-x2)^2 + (x2+x3-2)^2 + (x4-1)^2 + (x5-1)^2');
% Add Linear Constraints
B.AddCon('x1 + 3*x2 = 5');
B.AddCon('x3 + x4 - 2*x5 = 0');
B.AddCon('x2 - x5 = 0');

% Add Result Groups (the letter is arbitrary, any group name can be used)
B.AddResultGroup('A','Fuel Usage [l]');
B.AddResultGroup('B','Distance [km]');
% Add Result Expressions (Group:Name, Variable or Expression)
B.AddResultExp('A:Truck 1','x1');
B.AddResultExp('A:Truck 2','x2');
B.AddResultExp('B:Route 1','x3');
B.AddResultExp('B:Route 2','x4');
B.AddResultExp('B:Route 3','x5-x3'); %expressions with variables OK too

% Build Object
Build(B);

% Solve it and display formatted results
Solve(B,[ 2.5 0.5 2 -1 0.5 ]);
Results(B)

%% Example 6: Generating Code
% A new feature in OPTI v2.10 is the ability to generate C or C++ code from
% a nonlinear SymBuilder object. The default method is the generation of a
% MATLAB m-file with the nonlinear functions, however a C file may execute
% faster, and C++ files can use automatic derivatives via CppAD.

% C Code Example with Full (1st and 2nd) Analytical Derivatives
clc
% New SymBuilder Object
B = SymBuilder();
% Add Nonlinear Objective
B.AddObj('sin(pi*x1/12)*cos(pi*x2/16)');
% Add Linear Constraints
B.AddCon('-x1 + 2.5*x2 <= 1');
B.AddCon('x1 + 2.5*x2 <= -15');

% Build Object
Build(B);

% Solve Problem by Building a C callback function with analytical derivatives
sopts = symbset('cbmode','ccode');
[x,fval,ef,info] = Solve(B,[0;0],sopts)

% C++ Example with CppAD used to generate 1st and 2nd derivatives

% Solve Problem by Building a C++ callback compiled against CppAD
sopts = symbset('cbmode','cppad');
[x,fval,ef,info] = Solve(B,[0;0],sopts)

%% Example 7: Extracting Problem Data
% If you want to extract the functions or data generated by SymBuilder,
% simply use GetLinProb for LPs, GetQuadProb for QP/QCQPs and GetNLProb
% for NLPs to return an optiprob structure with the problem data.

clc
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

%% Example 8: SymBuilder Options
% SymBuilder uses the 'symbset' routine to control problem generation and
% solving. Use this routine to generate an options structure suitable for
% supplying to Solve or GetNLProb. 

clc
% Print all Options
symbset

% Set Solver, Display to 'iter'
sopts = symbset('solver','clp','display','iter');


%% Example 9: optisym
% A simple API to allow users familiar with fmincon type functions to use
% SymBuilder and leverage its functionality. Simply pass normal MATLAB
% functions to optisym and it internally converts them to Symbolic
% expressions and then to a SymBuilder object.
%
% optisym has the form:
%  [optiOpj,SymBobj] = optisym(fun,x0,lb,ub,con,cl,cu,xtype,sopts,verbose)
%
% As with normal SymBuilder models, if you pass a LP (even as function 
% handle), optisym will identify the problem as such and solve it as a LP. 
% This provides a convenient method to simplifying suitable models and
% solving them as efficiently as possible.
%
% Note this API is only compatible with functions that the Symbolic Toolbox
% can understand and parse.

% To illustrate, consider NLP1 Hock & Schittkowski #71
clc
% Normal Objective Function (min fun(x))
fun = @(x) x(1)*x(4)*sum(x(1:3)) + x(3);
% Nonlinear Constraint Function (cl <= con(x) <= cu)
con = @(x) [ prod(x);
             sum(x.^2)];
cl = [25;40];
cu = [Inf;40];
lb = ones(4,1);
ub = 5*ones(4,1);
x0 = [1 5 5 1]';
% Build OPTI Object
Opt = optisym(fun,x0,lb,ub,con,cl,cu)
% Solve
[x,fval,exitflag,info] = solve(Opt);
fval
info

% For more examples, see test_probs_sym.m in Test Problems/Development.




