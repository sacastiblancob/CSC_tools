# Compressed Sparse Column Tools 
> A series of tools developed to storage and manipulate matrices with CSC storage - MATLAB

## Installation

Just install matlab, and download the source codes, or by the next terminal command in a specific folder:

```sh
git clone https://github.com/sacastiblancob/CSC_tools .
```

## Usage example

To make some linear algebra you just need to use the CSC_tools.

If you want to have a regular matlab matriz A in CSC storage, you can use the function full2csc in the next way:

[Av,Ar,Ac] = full2csc(A)

But in fact that is not the correct way to do it, if you are making CFD, you are going to be able to come to your matrices through build them by diagonals and Kronecker products. Therefore, you can use tools as csc_kron or csc_diag to build matrices rather than convert regular matlab ones through full2csc. In fact, I strongly recommend to you to avoid the use of full2csc.

All the functions have a brief (sometimes detail) explanation of what that tool computes, what are the entries and the outputs. However, it is common to have a particular matrix, let's say A, stored in three different arrays:

 Av, Ar, Ac

Where Av have the values of the entries of A, Ar have the row indices of those entries, and Ac have the indices where every column start and the number of columns + 1 in the last position (because is CSC storage). And, for example, if you want to sum two matrices A+B=C (with csc_sum), it is intuitive to think that the entries will be Av, Ar, Ac, Bv, Br, Bc; and the outputs Cv, Cr, Cc. As it is. 

There are some functions that may not work as fast as I want to, but despite that, if you want to use csc_packlu to compute the LU decomposition of a matrix A, you should be aware of the structure of that matrix A and the structure of the result LU, because sometimes is more useful to have LU in full storage than in CSC storage. That depends on the application.

Be aware that for some solvers as Jacobi or SOR(Symmetric Over-Relaxation) you should use csc_prejacobi or csc_preSOR to prepare the entries for the actual solvers, which are csc_jacobi and csc_SOR. Those things are not together in the same function, because mostly you should solve a linear system of equations (sometimes several of them) every time step in a CFD simulation, so it is way better to compute csc_preSOR and csc_prejacobi before go into time loop.

## Folder Contents

./csc_CG --> Solves the system Ax=b with Conjugate Gradient.

./csc_diag --> Create a matrix in CSC by diagonals.

./csc_DPCG --> Solves the system Ax=b with Diagonal Preconditioned Conjugate Gradient method.

./csc_jacobi --> Solves the system Ax=b with Jacobi Iterative method (you should run csc_prejacobi before going into it)-

./csc_kron --> Computes the kronecker product between two matrices stored with this CSC metodology.

./csc_matmul --> Computes the matrix by matrix multiplication in CSC.

./csc_matvec --> Computes the matrix by vector multiplication in CSC.

./csc_prejacobi --> Prepares everything for Jacobi Iterative solver.

./csc_preSOR --> Prepares everything for Symmetric Over-Relaxation Iterative solver.

./csc_profs --> Script for testing.

./csc_SOR --> Solves the system Ax=b with Symmetric Over-Relaxation Iterative solver (you should run csc_preSOR before going into it).

./csc_spy --> Plots the sparsity of your matrix in CSC (the same as standar matlab spy()).

./csc_sum --> Computes the sum between two matrices in CSC.

./csc_trans --> Computes the transpose of a matrix in CSC.

./csc_packlu --> Computes the LU decomposition in packed way.

./csc_solpacklu --> Solves the system Ax=b with the result of LU decomposition.

./README.md --> You are standing here

## Meta
Sergio A. Castiblanco-Ballesteros

Bogota - Colombia

Mail 1: sacastiblancob@unal.edu.co

Mail 2: sergio.castiblanco@javeriana.edu.co

> Free Distribution and Open Source (As it should be!!)


