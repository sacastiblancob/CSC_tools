function [x,t,res] = csr_gmres(Av,Ac,Ar,b,x,m,niter,tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function solves the system Ax=b with the Restarted Generalized
% Minimal Residual (Restarted GMRES)
%
% Entries:
%     Av,Ac,Ar : Matrix A in CSR storage
%     b : Right hand side vector
%     x : First guest for the solution
%     m : Krylov space dimension
%     niter : Max. number of iterations
%     tol : tolerance for the stop through the norm of the residual
%
%      Sergio A. Castiblanco B.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = length(Ac)-1;
nb = length(b);
if n~=nb
    t = -1;
    return
end

for t = 1:niter

    %Residual
    r = b - csr_matvec(Av,Ac,Ar,x);
    
    %Arnoldi Method with Householder
    [HW,pbet] = csr_arnoldi_householder(Av,Ac,Ar,r,m);
    
    %Elimination process for Hessemberg matrix H through Givens Rotations
    g = zeros(m+1,1);
    g(1) = HW(1,1);
    for i=1:m
        d = sqrt(HW(i,i+1)^2 + HW(i+1,i+1)^2);
        s = HW(i+1,i+1)/d;
        c = HW(i,i+1)/d;
        HW(i,i+1) = HW(i,i+1)*c + HW(i+1,i+1)*s;
        HW(i+1,i+1) = 0;
        for j=i+1:m
            %Updating Hessemberg Matrix
            hij = HW(i,j+1);
            hi1j = HW(i+1,j+1);
            HW(i,j+1) = c*hij + s*hi1j;
            HW(i+1,j+1) = -s*hij + c*hi1j;
        end
        %Updating RHS vector g
        g(i+1) = -s*g(i);
        g(i) = c*g(i);
    end
    
    %Backward sustitution
    for i = m:-1:1
        g(i) = ( g(i)-HW(i,i+2:m+1)*g(i+1:m) )/HW(i,i+1);
    end
    %here y = g(1:m)
    
    %Computing z for xm = x0+z
    z = zeros(n,1);
    for k = m:-1:1
        z(k) = z(k) + g(k);
        q = z;
        bqTw = pbet(k)*(q(k:n)'*HW(k+1:n+1,k));
        for i = k:n
            z(i) = q(i) - HW(i+1,k)*bqTw;
        end
    end
    
    %updating the solution
    x = x + z;
    
    %tolerance over residual
    res = abs(g(m+1));
    if res<tol
        break
    end
end

end