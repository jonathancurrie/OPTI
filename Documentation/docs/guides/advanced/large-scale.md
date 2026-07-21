---
title: "Large-Scale Nonlinear Optimization"
slug: "/guides/advanced/large-scale/"
---

Solving large-scale (thousands of variables and constraints) nonlinear optimization problems do not require many changes in the way you pose the problem in MATLAB, but there are several techniques you can use to make solving them faster and more robust. This section will highlight a few tips and tricks for solving large-scale constrained nonlinear optimization problems ([NLPs](../../examples/problem-types/nlp.md) and [MINLPs](../../examples/problem-types/minlp.md)) using [IPOPT](../../solvers/ipopt.md) or [BONMIN](../../solvers/bonmin.md).

## General MATLAB Coding Guidelines
The techniques below will increase the speed of the solver, but are not designed increase the robustness of the problem description.

### Preallocating Memory
One of the simplest things you can do to improve performance of MATLAB code is to preallocate the memory required, before filling in array values. An example is shown below:

```matlab
%Inefficient Code
c(1,1) = 2*x(1);
c(2,1) = 3*x(2);
c(3,1) = 4*x(3);

%Efficent Code by Preallocating Memory
c = zeros(3,1);
c(1,1) = 2*x(1);
c(2,1) = 3*x(2);
c(3,1) = 4*x(3);
```

It is typical in optimization problems that you will know the size of your return matrix while writing the function (i.e. in design time), therefore preallocation is possible. If you do not preallocate your arrays MATLAB must fetch more memory everytime you create a new index in the array, then copy over the existing data. This process can be very time consuming for large problems.

NOTE: If you are using a white-box optimization solver such as SCIP or BARON, do not preallocate your return array! 

### Vectorizing Expressions
If you are a C programmer (like me) you may not always *think* in vectorized terms, but this is one of the most powerful features of MATLAB. Vectorized expressions almost always run faster so it is well worth the effort to start thinking in vectorized operations and functions. Consider the following example:

```matlab
%Inefficient (in MATLAB) C-like Code
j = 0;
for i = 1:length(x)
   j = j + x(i)*x(i);
end

%Efficent Code Using Vectorized Operators and Functions
j = sum(x.^2);
```

A general tip for vectorizing your code is to always ask at every loop: "Do I need this loop?" Typically `for` loops will be the easiest places to look for vectorized alternatives which may replace the loop altogether. Other areas include general algebraic equations. See if you can replace long equations with efficient linear algebra equivalents, such as for the common quadratic:

```matlab
%Explicitly Written Out Quadratic
j = 2*x(1)^2 + x(2)^2 + 4*x(1)*x(2) + 3*x(1);

%Same Quadratic using Matrix-Vector Expressions
H = [2 2; 2 1]; 
f = [3;0];
j = x'*H*x + f'*x;
```

In addition, using linear algebra often means your problem description scales easily as the problem size increases. Don't forget to use sparse matrices if you have less than around 10% non-zero entries.

### Use the Profiler
One of the most powerful tools for analyzing the performance of your functions (objective and constraints) is the MATLAB profiler. The profiler will run your program and automatically time and count all function and sub-function calls. The results are then displayed in an easy to navigate GUI which you can use to analyze the hot spots (computationally intensive regions of code). To initialize the profiler, use the following line of code:

```matlab
profile viewer
```

Once you have run your function from the profiler GUI, the results are automatically tabulated:

![profiler](/img/opti/profiler.png)

Using the results you can then manually optimize the hot spots within your functions, either using vectorization or rewriting sections of code. Profiling is one of the best tools for finding areas to increase the speed of your MATLAB code. Typically I obtain speed-ups in the range of 2-5x for large programs using the profiler and targeting hot spots with more efficient code.

### Building Large Sets of Nonlinear Constraints {#largenlcon}
A common problem people face is to how to create their nonlinear constraints, especially when they might have hundreds or thousands to define, or if the number changes frequently. If you have a set of equations that have the same structure, but vary in their indices between each equation, one way to automate this is to use a for-loop.

Consider the following series of quadratic equalities where *x* is the decision variable vector, and *r* is a constant vector:

![many nlcon](/img/opti/many_nlcon.png)

while it is possible to define this problem using matrix-vector notation, let's assume we don't want to. Obviously in the above problem, as N changes, the number of constraints change, as do the problem indices. Therefore one method to define this problem is to create a MATLAB function as follows

```matlab
function c = mynlcon(x,r)
%Nonlinear callback with adjustable number of constraints

N = length(r);

%Only preallocate if not solving with a white-box optimization solver
if(isa(x,'double'))
  c = zeros(N,1);
end

%Fill in each constraint equation (slow, but easy to read)
for i = 1:N
    c(i) = x(i)*x(2*N+i) - r(i)*x(i)*x(N+i) + (1-r(i))*x(N+i)*x(2*N+i);
end
```

Note that in this case it is possible to vectorize the expression, in which case, you should! This also allows us to create the function as a single anonymous function:

```matlab
%Vectorized version of the above
i = (1:N)';
nlcon = @(x) x(i).*x(2*N+i) - r(i).*x(i).*x(N+i) + (1-r(i)).*x(N+i).*x(2*N+i);
```

Both methods are equivalent when using a white-box optimization solver such as SCIP or BARON, however the vectorized version should be used with a general nonlinear solver (it should be faster). Whichever method you use, make sure you test your indexed expression first - write out the equations for a small N to ensure you are getting the correct answer!

## Supplying Derivatives
The techniques below may increase the speed of the solver, the robustness of the solver, or both.

### Provide Exact First Derivatives
When using a gradient based optimizer such as [IPOPT](../../solvers/ipopt.md) or [BONMIN](../../solvers/bonmin.md) relying on numerical approximations of the first derivatives, such as done by [mklJac](../../reference/utilities/djacobi.md) which uses centered finite-difference, may work acceptably for small, well posed problems, but it is unlikely to be robust enough for large, real-world problems. Not only will the problem be solved much slower when using numerical approximations, but it may not even converge if the gradient approximation is not accurate enough.

Providing exact (to numerical precision) first derivatives of the objective and constraints is typically the best thing you can do to help out a solver such as IPOPT in finding a solution. This documentation contains a page on [supplying first derivatives here](./deriv1.md).

### Provide Exact Second Derivatives
Supplying the Hessian of the Lagrangian (second derivatives of both the objective and all constraints) will save gradient based solvers such as [IPOPT](../../solvers/ipopt.md) or [BONMIN](../../solvers/bonmin.md) relying on quasi-newton numerical approximations of the second derivatives. For small problems supplying the Hessian may make no difference, but on large-scale problems I have seen it speed up the solver by factors up to 100x. The solver is also much more likely to converge to the (locally) optimal point with exact first and second derivatives.

This documentation contains a page on [supplying second derivatives here](./deriv2.md).

### Provide the Structure of the Derivative Matrices
Large-scale solvers typically require the structure of the Jacobian (first derivative of the constraints) and Hessian of the Lagrangian (second derivative of objective and constraints). The elements within these matrices may change value depending on **x**, thus simply supplying a function which calculates these elements does not tell the solver the structure (location of *all* non-zero elements) within these matrices. Based on this the solver must assume they are fully dense in order to account for every possible index having a value for every possible value of <small>**x**</small>, which is highly inefficient for large, sparse matrices as these typically are.

Therefore in order to leverage the inherent sparsity in these matrices you should always also supply the structure of these matrices via a sparse matrix of ones, where a one indicates a location of a non-zero element. Indicies not present within this structure matrix mean no matter what the value of <small>**x**</small> is, that index will never contain a value other than zero. This allows the solver to use efficient sparse linear algebra and linear solvers in order to speed up the internal computation of the solver. Note by default OPTI assumes these derivative matrices are fully dense, unless specified by the user.

See Example 5 in [supplying first derivatives here](./deriv1.md#expsparse) for more information on supplying the structure of these matrices.

### Ensure your Derivatives are Correct
When dealing with thousands of variables and constraints it is very easy to make a small mistake, perhaps have an index out by one or miss an element in one of the derivative matrices. Therefore in order to ensure your derivatives (first, second or both) are correct, OPTI provides a derivative checker which uses its own internal finite difference to check your supplied derivatives.

See Example 6 in [supplying first derivatives here](./deriv1.md#checkder) for more information on enabling OPTI's derivative checker. 

### Preallocating Sparse Derivative Matrices
As mentioned above, preallocating memory for arrays can save MATLAB repeatedly allocating memory as the array is built. The same applies for sparse matrices, which can be preallocated as follows:

```matlab
%Preallocate a 49x55 sparse matrix with 212 non-zero elements
J = spalloc(49,55,212);

%Fill in matrix
J(39,1) = 1;
J(40,1) = -1;
J(41,1) = x(1)/2.34;
J(42,1) = x(2)*x(3)*x(4);
J(31,2) = 1600;
J(48,2) = -1600;
%...
```

An alternative (and typically faster) method for building these matrices is to use sparse triplets, as described by Tim Davis in a MATLAB blog  [here](http://blogs.mathworks.com/loren/2007/03/01/creating-sparse-finite-element-matrices-in-matlab/). The article also describes how to create large, sparse, finite-element matrices, which is almost exactly what we are doing here! Well worth a read.

## Avoiding Derivatives
If your problem is (twice) differentiable, but entering all the derivative expressions manually does not appeal, there are alternative packages to use. Below are a few alternatives.

### Use Automatic Differentiation (AD)
Automatic differentiation is an alternative strategy to finite-difference (numerical approximation) and results in exact (to numerical precision) derivatives. OPTI provides a package called [ADIFF](../../reference/utilities/adiff.md) which is a simple implementation of AD in MATLAB for providing first derivatives. 

AD may take longer to evaluate than a finite-difference approximation, but the increased accuracy often results in less solver iterations, meaning you win in the long run. OPTI's ADIFF package only provides first derivatives for simple functions, so it is not quite ready yet for big problems where second derivatives are required.
 
### Use SymBuilder
SymBuilder is a small framework for creating algebraic optimization problems in MATLAB, and is included in the OPTI distribution. The advantage is that when the algebraic structure of the problem is retained, exact first and second derivatives can be automatically extracted using Symbolic techniques. This means you only have to enter your objective and constraint expressions, and SymBuilder will automatically generate all exact derivative functions for you! Note SymBuilder requires the MATLAB Symbolic Toolbox.

See the [SymBuilder documentation page](./sym-builder.md) for more information on using it.

### Use SCIP
The [SCIP](../../solvers/scip.md) interface released in OPTI v1.73 only requires the objective and constraint functions, as it internally generates the required first and second derivatives. This is done using  [CppAD](http://www.coin-or.org/CppAD/) , an automatic differentiation package in C++, which makes it very fast. However if you are not after a global solution, SCIP can take a long time to converge as it is a global solver.

### Use BARON
Support for [BARON](../../solvers/baron.md) was released in OPTI v2.05 and only requires the objective and constraint functions, as it internally generates the required first and second derivatives, much like SCIP does.

### Use a 3rd Party Package
Other than OPTI there are a number of alternative packages which aim to automate the generation of derivatives (among a lot of other functionality). I've listed a few below which you may like to try.

#### YALMIP [Free]
YALMIP is written and maintained by Johan Löfberg and allows high level model generation within MATLAB. Some of the solvers interfaced from YALMIP are supplied via OPTI, thus the packages work well together.

See the  [YALMIP documentation](https://yalmip.github.io/)  for more information.

#### BLOM [Free]
BLOM is a clever package which allows dynamic optimization within Simulink using IPOPT and other nonlinear solvers. 

See the  [BLOM Homepage](https://sites.google.com/berkeley.edu/mpc-lab/software/blom)  for more information.

#### AMPL [Commercial]
AMPL is a modelling language and presolver engine for a huge range of optimization problems. The resulting models can be solved either standalone from AMPL, or read using [amplRead](../../examples/file-formats/ampl.md) into MATLAB and solved using one of OPTI's solvers. When read using the AMPL interface into MATLAB, both first and second derivatives are also automatically added.

See the  [AMPL Homepage](http://www.ampl.com/)  for more information.

#### GAMS [Commercial]
GAMS is another modelling language and presolver engine for a huge range of optimization problems. GAMS includes a MATLAB interface but I have not had an opportunity to test it.

See the  [GAMS Homepage](https://www.gams.com/)  for more information.

#### TOMLAB [Commercial]
TOMLAB is a commercial package built on top of MATLAB, which links to many commercial solvers, in the much the same way OPTI links to open source solvers. However TOMLAB also includes an advanced automatic differentiation tool (MAD) which can provide exact first and second derivatives, including sparsity and structures.

See the  [TOMLAB Homepage](https://tomopt.com/tomlab/)  for more information.

### Use a Derivative-Free Solver
While some derivative free solvers work well for large scale problems, using a derivative-free solver simply because you can't be bothered entering the derivatives will typically always result in sub-optimal performance. If your problem is twice differentiable and convex, use IPOPT or BONMIN. As a rough rule of thumb only if your problem is not differentiable or non-convex should you consider a global / derivative-free solver.
