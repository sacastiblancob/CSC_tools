function [x,t,res] = csr_bicgstab(Av,Ac,Ar,b,x,r0a,niter,tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function solves the system Ax=b with the Biconjugate Gradient
% Stabilized Method (Van der Vorst algorithm).
%
% Entries:
%     Av,Ac,Ar : Matrix A in CSR storage
%     b : Right hand side vector
%     x : First guest for the solution
%     r0a : Arbitrary r0 for computing conjugate directions
%     niter : Max. number of iterations
%     tol : tolerance for the stop through the norm of the residual
%
%      Sergio A. Castiblanco B.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = b - csr_matvec(Av,Ac,Ar,x);
p = r;
for t=1:niter
    Ap = csr_matvec(Av,Ac,Ar,p);
    a = (r'*r0a)/(Ap'*r0a);
    s = r - a*Ap;
    As = csr_matvec(Av,Ac,Ar,s);
    w = (As'*s)/(As'*As);
    x = x + a*p + w*s;
    b = (a/w)*((s - w*As)'*r0a)/(r'*r0a);
    r = s - w*As;
    if norm(r)<tol
        break
    end
    p = r + b*(p - w*Ap);
end
res = norm(r);

end