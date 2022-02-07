function [x,t] = csr_CG(Av,Ac,Ar,b,x,niter,tol)
%
%This function solves the system Ax=b with the Conjugate Gradient method
%
% Entries:
%     Av,Ar,Ac : Matrix A in CSC storage
%     b : Right hand side vector
%     x : First guest for the solution
%     niter : Max. number of iterations
%     tol : tolerance for the stop through the norm of the residual
%
%
%      Sergio A. Castiblanco B. - Métodos Numéricos Avanzados
%      Pontificia Universidad Javeriana - Bogotá
%

ro = b - csr_matvec(Av,Ac,Ar,x);
d = ro;
for t = 1:niter
    Ad = csr_matvec(Av,Ac,Ar,d);
    if d'*Ad==0
        alf = 0.0;
    else
        alf = (ro'*ro)/(d'*Ad);
    end
    x = x + alf*d;
    if norm(b - csr_matvec(Av,Ac,Ar,x))<tol
        return
    end
    r = ro - alf*Ad;
    if ro'*ro==0
        bet = 0.0;
    else
        bet = (r'*r)/(ro'*ro);
    end
    d = r + bet*d;
    ro = r;
end

end