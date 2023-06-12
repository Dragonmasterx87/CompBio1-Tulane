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
mario <- c("mario", "luigi", "peaches", "donkeykong", "browser", "toad", 'giuseppe", "koopa", "spike", "penguinking", "kamek", "goomba")

# Check vector
mario

# Check if mario has any attributes assigned to it
attributes(mario)
```

## Names
Most common metadata assigned to vectors are names, dimensions and classes. You can look up each of these attributes using helper functions, which looks for that specific attribute.
```
# Check if mario has any name attributes assigned to it
names(mario)

# Lets assign names to mario vector
# Remember that the metadata wont change if you change the actual data.
names(mario) <- c("good", "good", "good", "good", "villan", "good", "good", villan", "villan", "good", "good", "villan") 

# Check if mario has any name attributes assigned to it
names(mario)

# R also displays metadata when you look at the object
mario
```

## Dims
You can transform an atomic vector into a complex n-dimensional array via the attribute `dim`. This can be especially useful when organising and classifying datasets.

```
# Lets check dim
dim(mario)

# Increase data dimensionality
dim(mario) <- (3, 4)

# Check data dimsnions now
dim(mario)
```
### Matrices
Matrices store data in a 2-dimensional array. Similar to dimensions in `dim` we can allocate the number of rows using the `nrow` command in the function `matrix`.

```
# Lets add 2 rows
mat <- matrix(mario, nrow = 2)
mat
```

### Array
This is a complex form of dataset organisation, where a dataset can be organised into a n-dimensional array. Think of a hypergeometric mode of data allocation, where you are storing data in the form of a complex hypercube of multiple dimensions. In this instance provide data in the first instance as a `vector` and dimensions as a `dim` in the the second instance.

```
# Lets create a hypergeometric array in the format of a hypercube
ar <- array(c(1:4, 5:8, 9:12, 13:16, ), dim = c(2, 2, 3, 3))
```

### Class
The class describes the format of data in R. You can check the `class` of data using the `class()` function.
```
class(mario)
class(ar)
```

If a dataset doesnt have a class assigned to it, the `class()` function will default to the vector type and output what type of vector data exists.
```
class("hello there!")
class(2019)
```

## Date & Time
We can check date and time that has passed from 12:00 a.m. Jan. 1, 1970 `POSIXct`
```
time_now <- Sys.time()
time_now

typeof(time_now)

class(time_now)
```

## Factors
Factors are R’s way of storing categorical information, like ethnicity or eye color. Think of a factor as something like a sex; it can only have certain values (male or female), and these values may have their own idiosyncratic order. This arrangement makes factors very useful for recording the treatment levels of a study and other categorical variables.

To make a factor, pass an atomic vector into the factor function. R will recode the data in the vector as integers and store the results in an integer vector. R will also add a levels attribute to the integer, which contains a set of labels for displaying the factor values, and a class attribute, which contains the class factor

gender <- factor(c("male", "female", "female", "male"))

```
typeof(gender)
## "integer"

attributes(gender)
## $levels
## [1] "female" "male"  
## 
## $class
## [1] "factor"
```















----



















