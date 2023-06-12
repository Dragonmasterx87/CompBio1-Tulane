---
title: Intro to R
layout: page
nav_order: 2
---

## Here we will go over a basic introduction to R

### What is R?
R is a free to download and use open-source computer language developed by Ross Ihaka and Robert Gentleman in 1993. It is an updated version of its predecessor the S programming language, which was commonly utilised in the early 1980s for statistical programming. R is designed to process and analyze complex statistical computations, store analysis and display results using complex visual graphics. R works in a sandbox known as the "environment"  which is a collection of objects such as variables, functions, tables etc. The underlying coding architecture of R is written in C, Fortran, and R itself. Complex workflows of functions are stored in the form of "packages" which can be installed in R and used to execute a series of tasks on a dataset. These packages are written in R however they may rely on C, C++, python and Fortran. R's computational philisophy is centered on object oriented programming (OOP).

### What is OOP and what is an "Object"?
In computer science an object can be a variable, data-structure, function, method or any data of a specific value, which can be classified by an identifier......so pretty much everything in R is an object. Objects can be of multiple different types, based on the nature of the variable. 

## Types of objects
### Atomic Vectors
This is the simplest form of vector data, where vector data is the simplest 1 dimensional array of any data type. There are 6 types of atomic vectors.

```
# Lets see different types of Atomic Vectors in R
# Atomic vector of type character.
print("abcdef")

# Atomic vector of type double.
print(89.214)

# Atomic vector of type integer.
print(93L)

# Atomic vector of type logical.
print(FALSE)

# Atomic vector of type complex.
print(9+36i)

# Atomic vector of type raw.
print(charToRaw('hi!'))
```
### Attributes
You can attach peices of info to atomic vectors, this can be information about the vectors itself and is helpful when classifying vector data. Think of this as "metadata"
```
# Lets make a vector
mario <- c("mario", "luigi", "peaches", "donkeykong", "browser", "toad", 'giuseppe", "koopa", "spike", "penguinking", "kamek", "goomba", "pirhanaplant")

# Check vector
mario

# Check if mario has any attributes assigned to it
attributes(mario)
```

## Names
Most common metadata assigned to vectors are names, dimensions and classes. You can look up each of these attributes using helper functions, which looks for that specific attribute.
```
# Check if mario has any name attributes assigned to it
names(mario) <- c("good", "good", "good", "good", "villan", "good", "good", villan", "villan", "good", "good", "villan", "villan") 

# Lets assign names to mario vector

```

----


