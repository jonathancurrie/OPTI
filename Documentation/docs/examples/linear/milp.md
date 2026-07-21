---
title: "Mixed Integer Linear Programming"
slug: "/examples/linear/milp/"
---

A Mixed Integer Linear Programming problem has the following form:

![eq MILP](/img/opti/eq_MILP.png)

Where **f** is a *n x 1* vector containing the linear objective function, which is subject to the following constraints: 

## Linear Inequalities*
**A** is a *m x n* sparse matrix, **b** is a *m x 1* vector 

## Linear Equalities*
**Aeq** is a *k x n* sparse matrix, **beq** is a *k x 1* vector 

## Decision Variable Bounds
**lb** and **ub** are *n x 1* vectors, where -inf or inf indicate an unbounded lower or upper bound, respectively 

## Integer Constraints
*x<sub>i</sub>* are decision variables which must be a integer number (...-2, -1, 0, 1, 2...) 

## Binary Constraints
*x<sub>j</sub>* are decision variables which must be a binary number (0,1), where *i ≠ j*. 

The goal is to minimize the objective function by selecting a value of **x** that also satisfies all constraints.
