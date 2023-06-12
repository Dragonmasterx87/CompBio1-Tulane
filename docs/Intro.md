---
title: Intro to R
layout: page
nav_order: 2
---

## Here we will go over a basic introduction to R

### What is R?
R is a free to download and use open-source computer language developed by Ross Ihaka and Robert Gentleman in 1993. It is an updated version of its predecessor the S programming language, which was commonly utilised in the early 1980s for statistical programming. R is designed to process and analyze complex statistical computations, store analysis and display results using complex visual graphics. R works in a sandbox known as the "environment"  which is a collection of objects such as variables, functions, tables etc. The underlying coding architecture of R is written in C, Fortran, and R itself. Complex workflows of functions are stored in the form of "packages" which can be installed in R and used to execute a series of tasks on a dataset. These packages are written in R however they may rely on C, C++, python and Fortran. R's computational philisophy is centered on object oreinted programming (OOP).



### What is OOP and what is an "Object"?
In computer science an object can be a variable, data-structure, function, method or any data of a specific value, which can be classified by an identifier......so pretty much everything in R is an object. Objects can be of multiple different types, based on the nature of the variable. 

## Types of objects
### Atomic Vectors
This is the simplest form of vector data, where vector data is the simplest 1 dimensional array of any data type. There are 6 types of atomic vectors.

```
# Atomic vector of type character.
print("abc")

# Atomic vector of type double.
print(12.5)

# Atomic vector of type integer.
print(63L)

# Atomic vector of type logical.
print(TRUE)

# Atomic vector of type complex.
print(2+3i)

# Atomic vector of type raw.
print(charToRaw('hello'))
```


----


