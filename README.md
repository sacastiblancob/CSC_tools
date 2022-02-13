# Compressed Sparse Column Tools / Now Also Compressed Sparse Row (and general CS)
> A series of tools developed to storage and manipulate matrices with CSC and CSR storage - MATLAB
Some are implemented in Fortran as well (I am working for the ones that are left).

## Installation

Just install matlab, and download the source codes, or by the next terminal command in a specific folder:

```sh
git clone https://github.com/sacastiblancob/CSC_tools .
```

## Usage example

To make some linear algebra you just need to use the CSS_tools (Compressed Sparse Storage Tools).

If you want to have a regular matlab matriz A in CSC storage or CSR storage, you can use the function full2csc in the next way:

[Av,Ar,Ac] = full2csc(A)
[Av,Ac,Ar] = full2csr(A)

But in fact that is not the correct way to do it, if you are making CFD, you are going to be able to come to your matrices through build them by diagonals and Kronecker products. Therefore, you can use tools as csc_kron or csc_diag to build matrices rather than convert regular matlab ones through full2csc. In fact, I strongly recommend to you to avoid the use of full2csc.

All the functions have a brief (sometimes detail) explanation of what that tool computes, what are the entries and the outputs. However, it is common to have a particular matrix, let's say A, stored in three different arrays:

 Av, Ar, Ac

Where Av have the values of the entries of A, Ar have the row indices of those entries, and Ac have the indices where every column start and the number of columns + 1 in the last position (because is CSC storage). And, for example, if you want to sum two matrices A+B=C (with csc_sum), it is intuitive to think that the entries will be Av, Ar, Ac, Bv, Br, Bc; and the outputs Cv, Cr, Cc. As it is. 

Or in CSR storage:

 Av, Ac, Ar

Where Av have the values of the entries of A, Ac have the column indices of those entries, and Ar have the indices where every row start and the number of columns + 1 in the last position (because is CSR storage). And, for example, if you want to sum two matrices A+B=C (with csr_symsum/csr_numsum), it is intuitive to think that the entries will be Av, Ac, Ar, Bv, Bc, Br; and the outputs Cv, Cc, Cr. As it is. 

There are some functions that may not work as fast as I want to, but despite that, if you want to use csc_packlu to compute the LU decomposition of a matrix A, you should be aware of the structure of that matrix A and the structure of the result LU, because sometimes is more useful to have LU in full storage than in CSC storage. That depends on the application. In fact, this algorithms (Direct Solvers) are not implemented yet in CSR storage, being in CSC badly-efficient.

Be aware that for some solvers as Jacobi or SOR(Symmetric Over-Relaxation) you should use csc_prejacobi or csc_preSOR to prepare the entries for the actual solvers, which are csc_jacobi and csc_SOR. Those things are not together in the same function, because mostly you should solve a linear system of equations (sometimes several of them) every time step in a CFD simulation, so it is way better to compute csc_preSOR and csc_prejacobi before go into time loop.

## Folder Contents

Both in ./fortran and ./matlab folders you will find almost the same functions/subroutines (listed below).
To use the Fortran subroutines you must copy the .f files in your fortran project source code folder, modify accordingly your Makefile (see examples in https://github.com/sacastiblancob/cavity-NS), and to call them, you just need to add the "USE CSC_STORAGE" statement in the declarations section of your fortran script.

./csc_CG ; csr_CG --> Solves the system Ax=b with Conjugate Gradient.

./csc_diag ; csr_diag --> Create a matrix in CSC by diagonals.

./csc_DPCG ; csr_DPCG --> Solves the system Ax=b with Diagonal Preconditioned Conjugate Gradient method.

./csc_jacobi ; csr_jacobi --> Solves the system Ax=b with Jacobi Iterative method (you should run csc_prejacobi before going into it)-

./css_kron --> Computes the kronecker product between two matrices stored with this Compressed Storage metodology (either CSC or CSR).

./csc_matmul --> Computes the matrix by matrix multiplication in CSC/CSR.

./csr_symmatmul ; csr_nummatmul --> Computes the matrix by matrix multiplication in CSR (symbolically and/or numerically).

./csc_matvec ; csr_matvec --> Computes the matrix by vector multiplication in CSC.

./csc_prejacobi ; csr_prejacobi --> Prepares everything for Jacobi Iterative solver.

./csc_preSOR ; csr_presor --> Prepares everything for Symmetric Over-Relaxation Iterative solver.

./csc_profs --> Script for testing.

./csc_SOR ; csr_SOR --> Solves the system Ax=b with Symmetric Over-Relaxation Iterative solver (you should run csc_preSOR before going into it).

./csc_spy --> Plots the sparsity of your matrix in CSC (the same as standar matlab spy()).

./csc_sum --> Computes the sum between two matrices in CSC (inneficient).

./csr_symsum ; csr_numsum --> Computes the sum between two matrices in CSR (symbolically and/or numerically).

./css_trans ; css_symtrans --> Computes the transpose of a matrix in Compressed Storage (symbolically and/or numerically).

./css_tperm --> Computes the transpose of matrix with permutations (apply it double for full permutations of row-columns).

./csc_packlu --> Computes the LU decomposition in packed way (inneficient).

./csc_solpacklu --> Solves the system Ax=b with the result of LU decomposition.

./csc_gmres ; csr_gmres --> Solves the system Ax=b with restarted Generalized Minimal Residual Method (restarted-GMRES).

./csc_bicgstab ; csr_bicgstab --> Solves the system Ax=b with Biconjugate Gradient Stabilized Method.

./housholderv --> Computes mirror vector for Householder Transformations (Golub & Van Loan algorithm).

./csc_arnoldi_householder ; csr_householder_arnoldi --> Computes Arnoldi's Method with Householder.

./csc_arnoldi_MGS ; csr_householder_MGS --> Computes Arnoldi's Method with Modified Gram-Schmidt.

./README.md --> You are standing here

## Meta
Sergio A. Castiblanco-Ballesteros

Bogota - Colombia

Mail 1: sacastiblancob@unal.edu.co

Mail 2: sergio.castiblanco@javeriana.edu.co

> Free Distribution and Open Source (As it should be!!)


