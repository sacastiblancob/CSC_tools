function [x,t] = csc_CG(Av,Ar,Ac,b,x,niter,tol,LUv,pc)
%
%This function solves the system Ax=b with the Conjugate Gradient method
%
% Entries:
%     Av,Ar,Ac : Matrix A in CSC storage
%     b : Right hand side vector
%     x : First guest for the solution
%     niter : Max. number of iterations
%     tol : tolerance for the stop through the norm of the residual
%     LUv  : values of LU decomposition matrix (for SSOR or ILU(0)) in
%              CSC_packed storage, if not SSOR or ILU(0), it must be have
%              arbitrary values.
%     pc  : preconditioning type
%       - 0 : No preconditioning
%       - 1 : Diagonal preconditioning
%       - 2 : Absolute diagonal preconditioning
%       - 3 : Symmetric SOR (Symmetric Gauss-Seidel -> SSOR, w=1)
%       - 4 : ILU(0)
%
%      Sergio A. Castiblanco B. - Métodos Numéricos Avanzados
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%No preconditioning
if pc==0

r = b - csc_matvec(Av,Ar,Ac,x);
p = r;
for t = 1:niter
    rr = r'*r;
    Ap = csc_matvec(Av,Ar,Ac,p);
    a = rr/(Ap'*p);
    x = x + a*p;
    r = r - a*Ap;
    if norm(r)<tol
        break
    end
    b = (r'*r)/rr;
    p = r + b*p;
end

%Diagonal preconditioning
elseif pc==1

r = b - csc_matvec(Av,Ar,Ac,x);
D = csc_diaga(Av,Ar,Ac);
z = r./D;
p = z;
for t=1:niter
    Ap = csc_matvec(Av,Ar,Ac,p);
    rz = r'*z;
    a = (rz)/(Ap'*p);
    x = x + a*p;
    r = r - a*Ap;
    if norm(r)<tol
        break
    end
    z = r./D;
    b = (r'*z)/(rz);
    p = z + b*p;
end

%Absolute Diagonal preconditioning
elseif pc==2

r = b - csc_matvec(Av,Ar,Ac,x);
D = abs(csc_diaga(Av,Ar,Ac));
z = r./D;
p = z;
for t=1:niter
    Ap = csc_matvec(Av,Ar,Ac,p);
    rz = r'*z;
    a = (rz)/(Ap'*p);
    x = x + a*p;
    r = r - a*Ap;
    if norm(r)<tol
        break
    end
    z = r./D;
    b = (r'*z)/(rz);
    p = z + b*p;
end

%SSOR or ILU preconditioning
elseif pc==3 || pc==4

r = b - csc_matvec(Av,Ar,Ac,x);
z = csc_solpacklu(LUv, Ar, Ac, r);
p = z;
for t=1:niter
    Ap = csc_matvec(Av,Ar,Ac,p);
    rz = r'*z;
    a = (rz)/(Ap'*p);
    x = x + a*p;
    r = r - a*Ap;
    if norm(r)<tol
        break
    end
    z = csc_solpacklu(LUv, Ar, Ac, r);
    b = (r'*z)/(rz);
    p = z + b*p;
end

end

end