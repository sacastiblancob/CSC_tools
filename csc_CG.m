function [x,t] = csc_CG(Av,Ar,Ac,b,x,niter,tol)
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

ro = b - csc_matvec(Av,Ar,Ac,x);
d = ro;
for t = 1:niter
    Ad = csc_matvec(Av,Ar,Ac,d);
    alf = (ro'*ro)/(d'*Ad);
    x = x + alf*d;
    if norm(b - csc_matvec(Av,Ar,Ac,x))<tol
        return
    end
    r = ro - alf*Ad;
    bet = (r'*r)/(ro'*ro);
    d = r + bet*d;
    ro = r;
end

end