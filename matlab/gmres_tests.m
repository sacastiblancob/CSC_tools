% %
% % Full Orthogonalization Method (FOM) and related things
% %
% % Full Orthogonalization Method (FOM) with Modified Gram-Schmidt in CSC
% h = [1 3 1;1 3 1;2 6 2;2 6 2;3 8 3;3 8 3];
% A = diag(h(1:5,1),-1) + diag(h(:,2),0) + diag(h(1:5,3),1);
% [Av,Ar,Ac] = full2csc(A);
% b = ones(length(A),1);
% x0 = ones(length(A),1);
% xe = A\b;
% ni = n;
% % E = zeros(ni,1);
% 
% % for i=1:ni
%     r0 = b - csc_matvec(Av,Ar,Ac,x0);
%     bet = norm(r0);
%     v = r0/bet;
%     m = 6;
%     [V,H] = arnoldi_MGS_csc(Av,Ar,Ac,v,m,eps(1E6));
%     bv = zeros(m,1);
%     bv(1) = bet;
%     y = inv(H)*bv;   %This system H*y=bv should be computed with Givens Rotations instead
%     x0 = x0 + V*y;
% %     E(i) = norm(x0-xe);
% % end
% % semilogy(E)
% % As it becomes expensive when m -> n, the restarted version is proposed.
% % In which m is set fixed or expanded up to a m_max

% % Restarded-FOM with Modified Gram-Schmidt in CSC
% h = [1 3 1;1 3 1;2 6 2;2 6 2;3 8 3;3 8 3];
% A = diag(h(1:5,1),-1) + diag(h(:,2),0) + diag(h(1:5,3),1);
% [Av,Ar,Ac] = full2csc(A);
% n = length(A);
% b = ones(length(A),1);
% x0 = ones(length(A),1);
% xe = A\b;
% m = 3;
% ni = n;
% E = zeros(ni,1);
% 
% for i=1:ni
%     r0 = b - csc_matvec(Av,Ar,Ac,x0);
%     bet = norm(r0);
%     v = r0/bet;
%     [V,H] = arnoldi_MGS_csc(Av,Ar,Ac,v,m,eps(1E6));
%     bv = zeros(m,1);
%     bv(1) = bet;
%     y = inv(H)*bv;   %This system H*y=bv should be computed with Givens Rotations instead
%     x0 = x0 + V*y;
%     E(i) = norm(x0-xe);
% end
% semilogy(E2)
% hold on
% semilogy(E3)
%
% %
% % Generalized Minimal Residual (GMRES) and related things
% %
% 
% %generic case
% h = [1 3 1;1 3 1;2 6 2;2 6 2;3 8 3;3 8 3];
% A = diag(h(1:5,1),-1) + diag(h(:,2),0) + diag(h(1:5,3),1);
% % A(1,5) = -10;
% % A(6,2) = 5;
% [Av,Ar,Ac] = full2csc(A);
% n = length(A);
% b = ones(n,1);
% xe = A\b;
% 
% % %Finite differences 1d-Poisson Equation with neuman boundary conditions
% % % i.e, A is singular
% % %
% % % Here, the problem, after regularization is solved by GMRES exactly if one
% % % uses a krylov space of dimmension m = n/2
% % %
% % n = 10;
% % A1 = -ones(n-1,1);
% % A0 = 2*ones(n,1);
% % A0(1) = 1;
% % A0(n) = 1;
% % % h = [-1 1 -1;-1 2 -1;-1 2 -1;-1 2 -1;-1 2 -1;-1 1 -1];
% % A = diag(A1,-1) + diag(A0,0) + diag(A1,1);
% % [Av,Ar,Ac] = full2csc(A);
% % n = length(A);
% % b = ones(n,1);
% % b(1) = 0;
% % b(n) = 0;
% % xe = A\b;
% 
% % %regularization
% % [UL,S,V] = svds(A,1,'smallest');
% % RM = (eye(n)-UL*UL');
% % b = RM*b;
% 
% % %arbitrary x0
% x0 = ones(n,1);
% r0 = b - csc_matvec(Av,Ar,Ac,x0);
% 
% % %What if r0 = eigvector %nays!! in this case with m=1 is enough
% % [V,D] = eigs(A,2);
% % r0 = 1*V(:,1) + 2*V(:,2);
% % x0 = b-r0;
% % x0 = A\x0;
% % % Breakedowns emerges when r_0 is a lineal combination of span{(eigv_A)_k}
% % % (any k eigenvectors of A), if that is the case, when m=k, the solution of
% % % Ax=b is exact. And breakedowns emerge for steps m=>j>k. This is clear for
% % % real SPD matrices at least.
% 
% %Arnoldi Method with Householder
% m = 5;
% [HW,pbet] = csc_arnoldi_householder(Av,Ar,Ac,r0,m);
% % HWp = HW;
% 
% %
% %Solving H*y = (beta*e_1) for y
% 
% %Elimination process for Hessemberg matrix H through Givens Rotations
% % i = 1;
% g = zeros(m+1,1);
% g(1) = HW(1,1);
% for i=1:m
%     d = sqrt(HW(i,i+1)^2 + HW(i+1,i+1)^2);
%     s = HW(i+1,i+1)/d;
%     c = HW(i,i+1)/d;
%     HW(i,i+1) = HW(i,i+1)*c + HW(i+1,i+1)*s;
%     HW(i+1,i+1) = 0;
%     for j=i+1:m
%         %Updating Hessemberg Matrix
%         hij = HW(i,j+1);
%         hi1j = HW(i+1,j+1);
%         HW(i,j+1) = c*hij + s*hi1j;
%         HW(i+1,j+1) = -s*hij + c*hi1j;
%     end
%     %Updating RHS vector g
%     g(i+1) = -s*g(i);
%     g(i) = c*g(i);
% end
% % ge = g;
% 
% %Backward sustitution
% for i = m:-1:1
%     g(i) = ( g(i)-HW(i,i+2:m+1)*g(i+1:m) )/HW(i,i+1);
% end
% %here y = g(1:m)
% 
% %Computing z for xm = x0+z
% z = zeros(n,1);
% for k = m:-1:1
%     z(k) = z(k) + g(k);
%     q = z;
%     bqTw = pbet(k)*(q(k:n)'*HW(k+1:n+1,k));
%     for i = k:n
%         z(i) = q(i) - HW(i+1,k)*bqTw;
%     end
% end
% 
% x1 = x0 + z;
% 
% %this computation is not needed since r1 = abs(g(m+1))
% r1 = norm(b - csc_matvec(Av,Ar,Ac,x1));

% %
% % Generalized Minimal Residual (GMRES) Function
% %

% %generic case
% h = [1 4 1;1 7 1;2 8 2;2 8 2;3 10 3;3 10 3];
% A = diag(h(1:5,1),-1) + diag(h(:,2),0) + diag(h(1:5,3),1);
% A(1,5) = -7;
% A(5,1) = 7;
% [Av,Ar,Ac] = full2csc(A);
% n = length(A);
% b = ones(n,1);
% xe = A\b;
% 
% %solving with GMRES
% m = 2;
% niter = 100;
% tol = 1E-10;
% x = b;
% [xgr,tgr,res] = csc_gmres(Av,Ar,Ac,b,x,m,niter,tol);
% [xcg,tcg] = csc_CG(Av,Ar,Ac,b,x,niter,tol);
% w=0.9;
% [Pv,Pr,Pc,Qv,Qr,Qc] = csc_preSOR(Av,Ac,Ar,w);
% [xsor,tsor] = csc_SOR(Av,Ar,Ac,Pv,Pr,Pc,Qv,Qr,Qc,b,x,niter,tol);

% A = [3 -4.2;1 3];
% b = A*[1;1];
% [Av,Ar,Ac] = full2csc(A);
% x = [0;-1];
% m=2;
% niter=50;
% tol = 1E-10;
% [x,t,res] = csc_gmres(Av,Ar,Ac,b,x,m,niter,tol);

%
% Biconjugate Gradient Stabilized Method (BICGSTAB)
%
h = [1 4 1;1 7 1;2 8 2;2 8 2;3 10 3;3 10 3];
A = diag(h(1:5,1),-1) + diag(h(:,2),0) + diag(h(1:5,3),1);
A(3,6) = 30;
A(3,1) = 20;
[Av,Ar,Ac] = full2csc(A);
n = length(A);
% b = ones(n,1);
b = [1;2;3;4;5;6];
xe = A\b;
niter = 5;
tol = 1E-10;

x = b;
m = 6;
[xg,tg,rg] = csc_gmres(Av,Ar,Ac,b,x,m,niter,tol);

r = b - csc_matvec(Av,Ar,Ac,x);
% r0a = r;
% r0a = A'*r;
% r0a = A'*(b - A'*r);
r0a = csc_matvec(Av,Ar,Ac,r); % Among the options related to first r, this 
                              % gives the best performance for now for r0 arbitrary
ep = zeros(n,1);
ep(4) = 1;
% ep = ones(n,1);
%the unitary vector ep(p) works fine for r0 arbitrary, but is not clear for which p
%depending on matrix A, however best p does not seem to depends on vector b
[xb,tb,rb] = csc_bicgstab(Av,Ar,Ac,b,x,ep,niter,tol);

% r0a = r;
% p = r;
% for t=1:niter
%     Ap = csc_matvec(Av,Ar,Ac,p);
%     a = (r'*r0a)/(Ap'*r0a);
%     s = r - a*Ap;
%     As = csc_matvec(Av,Ar,Ac,s);
%     w = (As'*s)/(As'*As);
%     x = x + a*p + w*s;
%     b = (a/w)*((s - w*As)'*r0a)/(r'*r0a);
%     r = s - w*As;
%     if norm(r)<tol
%         break
%     end
%     p = r + b*(p - w*Ap);
% end
















