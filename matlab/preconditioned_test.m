%
% Preconditioned Solvers
%
% CG, BICGSTAB, GMRES
% -> Diagonal
% -> Absolute Diagonal
% -> SOR/SSOR
% -> ILU(0)
% -> ILU(1) ??

% %
% %CONJUGATE GRADIENT
% %
% %generic case
% % h = [3 6 3;3 6 3;6 8 6;6 8 6;8 10 8;8 10 8];
% % h = [3 -6 3;3 -6 3;6 8 6;6 -8 6;8 10 8;8 -10 8];
% % A = diag(h(1:5,1),-1) + diag(h(:,2),0) + diag(h(1:5,3),1);
% h = [3 -6 3;3 -6 3;6 8 6;6 -8 6;8 10 8;8 -10 8];
% A = diag(h(1:4,1),-2) + diag(h(:,2),0) + diag(h(1:4,3),2);
% A(5,3) = 0;
% A(3,5) = 0;
% mE = tril(A) - diag(diag(A));
% Da = diag(diag(A));
% mF = triu(A) - diag(diag(A));
% [Av,Ar,Ac] = full2csc(A);
% [Bv,Bc,Br] = full2csr(A);
% LUv = zeros(length(Av));
% n = length(A);
% b = ones(n,1);
% xe = A\b;
% niter = 6;
% tol = 1E-10;
% 
% %No preconditioning
% [xc,tc] = csc_CG(Av,Ar,Ac,b,b,niter,tol,LUv,0);
% [xcr,tcr] = csr_CG(Bv,Bc,Br,b,b,niter,tol,LUv,0);
% nc = norm(xc - xe);
% ncr = norm(xcr - xe);
% 
% %Diagonal preconditioned CG
% [xcd,tcd] = csc_CG(Av,Ar,Ac,b,b,niter,tol,LUv,1);
% [xcdr,tcdr] = csr_CG(Bv,Bc,Br,b,b,niter,tol,LUv,1);
% ncd = norm(xcd - xe);
% ncdr = norm(xcdr - xe);
% 
% 
% %Absolute Diagonal preconditioned CG
% [xcad,tcad] = csc_CG(Av,Ar,Ac,b,b,niter,tol,LUv,2);
% [xcadr,tcadr] = csr_CG(Bv,Bc,Br,b,b,niter,tol,LUv,2);
% ncad = norm(xcad - xe);
% ncadr = norm(xcadr - xe);
% 
% %Symmetric SOR
% w = 1;
% [PQv] = csc_preconSSOR(Av,Ar,Ac,w);
% [PQvr] = csr_preconSSOR(Bv,Bc,Br,w);
% [xcs,tcs] = csc_CG(Av,Ar,Ac,b,b,niter,tol,PQv,3);
% [xcsr,tcsr] = csr_CG(Bv,Bc,Br,b,b,niter,tol,PQvr,3);
% ncs = norm(xcs - xe);
% ncsr = norm(xcsr - xe);
% 
% %ILU(0)
% [LUv] = csc_preconILU0(Av,Ar,Ac);
% [LUvr] = css_trans(LUv,Ar,Ac);
% [xcl,tcl] = csc_CG(Av,Ar,Ac,b,b,niter,tol,LUv,4);
% [xclr,tclr] = csr_CG(Bv,Bc,Br,b,b,niter,tol,LUvr,4);
% ncl = norm(xcl - xe);
% nclr = norm(xclr - xe);
% 
% % % x0 = b;
% % % r = b - csc_matvec(Av,Ar,Ac,x0);
% % % P = eye(n) + w*mE*inv(Da);
% % % Q = Da + w*mF;
% % % PQ = (tril(P) - diag(diag(P))) + Q;
% % % yr = inv(P)*r;
% % % zr = inv(Q)*yr;
% % % z = csc_solpacklu(PQv, Ar, Ac, r);
% % % z2 = csr_solpacklu(PQvr,Acr,Arr, r);

% %
% %BICONJUGATE GRADIENT
% %
% %generic case
% % h = [3 6 3;3 6 3;6 8 6;6 8 6;8 10 8;8 10 8];
% h = [3 6 -3;4 6 3;-10 8 6;6 8 8;-8 -10 -8;8 -1 4];
% A = diag(h(1:5,1),-1) + diag(h(:,2),0) + diag(h(1:5,3),1);
% A(1,5) = -4;
% A(6,2) = -8;
% A(6,1) = 15;
% [Av,Ar,Ac] = full2csc(A);
% [Bv,Bc,Br] = full2csr(A);
% LUv = zeros(length(Av));
% n = length(A);
% b = ones(n,1);
% xe = A\b;
% niter = n-1;
% tol = 1E-10;
% 
% % No preconditioning
% x0 = b;
% r0a = b - csc_matvec(Av,Ar,Ac,x0);
% [xbi,tbi,resbi] = csc_bicgstab(Av,Ar,Ac,b,b,r0a,niter,tol,LUv,0);
% [xbir,tbir,resbir] = csr_bicgstab(Bv,Bc,Br,b,b,r0a,niter,tol,LUv,0);
% nbi = norm(xbi - xe);
% nbir = norm(xbir - xe);
% 
% % Diagonal Preconditioning
% [xbid,tbid,resbid] = csc_bicgstab(Av,Ar,Ac,b,b,r0a,niter,tol,LUv,1);
% [xbidr,tbidr,resbidr] = csr_bicgstab(Bv,Bc,Br,b,b,r0a,niter,tol,LUv,1);
% nbid = norm(xbid - xe);
% nbidr = norm(xbidr - xe);
% 
% % Absolute Diagonal Preconditioning
% [xbiad,tbiad,resbiad] = csc_bicgstab(Av,Ar,Ac,b,b,r0a,niter,tol,LUv,2);
% [xbiadr,tbiadr,resbiadr] = csr_bicgstab(Bv,Bc,Br,b,b,r0a,niter,tol,LUv,2);
% nbiad = norm(xbiad - xe);
% nbiadr = norm(xbiadr - xe);
% 
% % SSOR Preconditioning
% w = 0.9;
% [PQv] = csc_preconSSOR(Av,Ar,Ac,w);
% [PQvr] = csr_preconSSOR(Bv,Bc,Br,w);
% [xbis,tbis,resbis] = csc_bicgstab(Av,Ar,Ac,b,b,r0a,niter,tol,PQv,3);
% [xbisr,tbisr,resbisr] = csr_bicgstab(Bv,Bc,Br,b,b,r0a,niter,tol,PQvr,3);
% nbis = norm(xbis - xe);
% nbisr = norm(xbisr - xe);
% 
% % ILU(0) Preconditioning
% [LUv] = csc_preconILU0(Av,Ar,Ac);
% [LUvr] = css_trans(LUv,Ar,Ac);
% [xbil,tbil,resbil] = csc_bicgstab(Av,Ar,Ac,b,b,r0a,niter,tol,LUv,4);
% [xbilr,tbilr,resbilr] = csr_bicgstab(Bv,Bc,Br,b,b,r0a,niter,tol,LUvr,4);
% nbil = norm(xbil - xe);
% nbilr = norm(xbilr - xe);

% 
% GMRES
%
% h = [3 6 3;3 6 3;6 8 6;6 8 6;8 10 8;8 10 8];
% A = diag(h(1:5,1),-1) + diag(h(:,2),0) + diag(h(1:5,3),1);
h = [3 6 -3;4 6 3;-10 8 6;6 8 8;-8 -10 -8;8 -1 4];
A = diag(h(1:5,1),-1) + diag(h(:,2),0) + diag(h(1:5,3),1);
A(1,5) = -4;
A(6,2) = -8;
A(6,1) = 15;
[Av,Ar,Ac] = full2csc(A);
[Bv,Bc,Br] = full2csr(A);
LUv = zeros(length(Av));
n = length(A);
b = ones(n,1);
xe = A\b;
niter = 100;
tol = 1E-10;
m = n-5;
x = b;

% No Preconditioning
[xg,tg,resg] = csc_gmres(Av,Ar,Ac,b,b,m,niter,tol,LUv,0);
[xgr,tgr,resgr] = csr_gmres(Bv,Bc,Br,b,b,m,niter,tol,LUv,0);
ng = norm(xg-xe);
ngr = norm(xgr-xe);

% Diagonal Preconditioning
[xgd,tgd,resgd] = csc_gmres(Av,Ar,Ac,b,b,m,niter,tol,LUv,1);
[xgdr,tgdr,resgdr] = csr_gmres(Bv,Bc,Br,b,b,m,niter,tol,LUv,1);
ngd = norm(xgd-xe);
ngdr = norm(xgdr-xe);

% Absolute Diagonal Preconditioning
[xgad,tgad,resgad] = csc_gmres(Av,Ar,Ac,b,b,m,niter,tol,LUv,2);
[xgadr,tgadr,resgadr] = csr_gmres(Bv,Bc,Br,b,b,m,niter,tol,LUv,2);
ngad = norm(xgad-xe);
ngadr = norm(xgadr-xe);

% SSOR Preconditioning
w = 1.1;
[PQv] = csc_preconSSOR(Av,Ar,Ac,w);
[PQvr] = csr_preconSSOR(Bv,Bc,Br,w);
[xgs,tgs,resgs] = csc_gmres(Av,Ar,Ac,b,b,m,niter,tol,PQv,3);
[xgsr,tgsr,resgsr] = csr_gmres(Bv,Bc,Br,b,b,m,niter,tol,PQvr,3);
ngs = norm(xgs-xe);
ngsr = norm(xgsr-xe);

% ILU(0) Preconditioning
[LUv] = csc_preconILU0(Av,Ar,Ac);
[LUvr] = css_trans(LUv,Ar,Ac);
[xgl,tgl,resgl] = csc_gmres(Av,Ar,Ac,b,b,m,niter,tol,LUv,4);
[xglr,tglr,resglr] = csr_gmres(Bv,Bc,Br,b,b,m,niter,tol,LUvr,4);
ngl = norm(xgl-xe);
nglr = norm(xglr-xe);

% % ILU(0)
% h = [3 -6 -2;3 6 -2;6 -8 -4;6 8 -4;8 -10 -6;8 10 -6];
% A = diag(h(1:5,1),-1) + diag(h(:,2),0) + diag(h(1:5,3),1);
% A(1,5) = -4;
% A(6,2) = -8;
% [Av,Ar,Ac] = full2csc(A);
% % [Bv,Bc,Br] = full2csr(A);
% % [Av,Ac,Ar] = full2csr(A);
% LUv = zeros(length(Av));
% n = length(A);
% b = ones(n,1);
% xe = A\b;
% niter = 10;
% tol = 1E-10;
% m = 4;
% x = b;
% 
% % %ILU(0) naive
% % % A2 = diag(diag(A)) + A;
% % % ILU0 = A2;
% % ILU0 = A;
% % 
% % %packed ILU0 re Naive
% % for i=2:n
% %     for k=1:i-1
% %         if (A(i,k)~=0)
% %             ILU0(i,k) = ILU0(i,k)/ILU0(k,k);
% %             for j=k+1:n
% %                 if (A(i,j)~=0)
% %                     ILU0(i,j) = ILU0(i,j) - ILU0(i,k)*ILU0(k,j);
% %                 end
% %             end
% %         end
% %     end
% % end
% % 
% % %retrieving L and U
% % L0 = eye(n);
% % U0 = zeros(n,n);
% % for i=1:n
% %     for j=1:n
% %         if i<=j
% %             U0(i,j) = ILU0(i,j);
% %         else
% %             L0(i,j) = ILU0(i,j);
% %         end
% %     end
% % end
% % R2 = L0*U0 - A;
% % 
% % %init L and U, y and x
% % L = eye(n);
% % U = A;
% % Aorg = A;
% % y = zeros(n,1);
% % x = zeros(n,1);
% % 
% % for j = 1:n-1
% %     % Subtract multiples lj of A(j,:) from rows that follow
% %     lj = U(j+1:n,j)/U(j,j);
% %     L(j+1:n,j) = lj;
% %     U(j+1:n,:) = U(j+1:n,:)-lj*U(j,:);
% % end
% 
% %Incomplete LU sparse CSC
% % LUv = Av;
% % m = length(Ac)-1;
% % nz = length(Av);
% % for j=1:m-1
% %     for i=Ac(j):Ac(j+1)-1
% %         if j == Ar(i)
% %             d = LUv(i);
% %         end
% %         if j < Ar(i)
% %             LUv(i) = LUv(i)/d;
% %             for j2 = j+1:m
% %                 a = 0;
% %                 for i2 = Ac(j2):Ac(j2+1)-1
% %                     if Ar(i2)==j
% %                        a = LUv(i2);
% %                     end
% %                     if Ar(i)==Ar(i2)
% %                         LUv(i2) = LUv(i2) - LUv(i)*a;
% %                     end
% %                 end
% %             end            
% %         end
% %     end
% % end
% [LUv] = csc_preconILU0(Av,Ar,Ac);
% [LUvr,Acr,Arr] = css_trans(LUv,Ar,Ac);







