% This code shows didactically One-dimensional Projection Iterative Methdos
% for solving linear system of equations.
% Namely:
%   Steepest Descendt
%   Minimal Residual Iteration
%   Residual Norm Steepest Descent
%
% Some theory:
%
% For the problem:
%               Ax = b
% Where A is a non-singular square matrix of dimensions nxn.
%
% Projection Iterative Methods works by choosing two spaces K and L, where:
%   K -> Subspace of candidate approximants, or search subspace of
%   dimension m < n
%   L -> Subspace of constraints or left subspace of dimension m
%
%   Then,
%       - m restrictions must be settled.
%       - Tipical way: m orthogonal conditions.
%       - Specifically: Residual vector b-Ax is constrained to be
%       orthogonal to L subspace.
%
% When K = L, Galerkin case (Orthogonal)
% When K /= L, Petrov-Galerkin case (Oblique)
%
% Conditions
%       x_1 = x_0 + d   ,   d € K               (1)
%     (r_o - Ad, w)=0   ,   for all w € L       (2)
%
% Example: Gauss-seidel can be seen as an projection method, in which every
% step, K=L=span{e_i}, for i cycling from i=1 to i=1,...,n
%
% One-Dimensional cases:
%       K = span{v}   ;   L = span{w}
% Where v and w are two vectors.
%
% From condition (1):
%   x_1 = x_0 + d
% As, d € K, then, d must be a scalar product of vector v. Then:
%   x_1 = x_0 + a*v
% By condition (2), a must be equal to:
%   a = (r,w)/(Av,w)
% Where the operator (,) is the inner-product. Now, how to choose v and w?
%
% Steepest descend --> v=r   ;   w=r  (Galerki condition using v=r)
% Minimal Residual Iteration --> v=r   ;   w=Ar (Petrov-Galerkin using v=r)
% Residual Norm Steepest Descendt --> v=A'r   ;   w=Av (Petrov-Galerkin using v=A'r)
%
% Conditions over A:
% Steepest Descend -> A must be non-singular Symmetric Positive Definite.
% Minimal Residual Iteration -> A must be non-singular Positive Definite.
% Residual Norm Steepest Descendt -> A must be non-singular.
%
% For further details and other additional theoretical results, see:
%   Yousef Saad - Iterative Methods for Sparse Linear Systems
%
% Made by Sergio C.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Steepest Descend (s) ; Minimal Residual (m) ; Residual Norm SD (r)

%Matrix A, b, and exact solution
A = [3 -4.2;1 3];
b = A*[1;1];
xe = A\b;  %exact solution

%Initial guess
x = [0;-1];

%Algorithms
ni = 10;    %Number of iterations

r = b - A*x;
p = A*r;

rs = r; ps = p; xs = x;        %steepest descend
rm = r; pm = p; xm = x;        %Minimal Residual
rr = r; pr = p; xr = x;        %Residual Norm Steepest Descend

%For storing process
%steepest descend
Rs = zeros(length(b),ni+1);
Ps = zeros(length(b),ni+1);
Xs = zeros(length(b),ni+1);
Es = zeros(length(b),ni+1);
As = zeros(1,ni+1);
Rs(:,1) = rs;
Ps(:,1) = ps;
Xs(:,1) = xs;
Es(:,1) = xe - xs;

%Minimal Residual
Rm = zeros(length(b),ni+1);
Pm = zeros(length(b),ni+1);
Xm = zeros(length(b),ni+1);
Em = zeros(length(b),ni+1);
Am = zeros(1,ni+1);
Rm(:,1) = rm;
Pm(:,1) = pm;
Xm(:,1) = xm;
Em(:,1) = xe - xm;

%Residual Norm Steepest Descend
Rr = zeros(length(b),ni+1);
Vr = zeros(length(b),ni+1);
Wr = zeros(length(b),ni+1);
Xr = zeros(length(b),ni+1);
Er = zeros(length(b),ni+1);
Ar = zeros(1,ni+1);
Rr(:,1) = rr;
% Vr(:,1) = A'*rr;
% Wr(:,1) = A*Vr(:,1);
Xr(:,1) = xr;
Er(:,1) = xe - xr;

for t=1:ni
    %
    %Steepest Descend
    %
    as = (rs'*rs)/(ps'*rs);
    xs = xs + as*rs;
    rs = rs - as*ps;
    ps = A*rs;
    
    %updating
    Rs(:,t+1) = rs;
    Ps(:,t+1) = ps;
    Xs(:,t+1) = xs;
    Es(:,t+1) = xe - xs;
    As(:,t) = as;

    %
    %Minimal Residual Iteration
    %
    am = (pm'*rm)/(pm'*pm);
    xm = xm + am*rm;
    rm = rm - am*pm;
    pm = A*rm;

    %updating
    Rm(:,t+1) = rm;
    Pm(:,t+1) = pm;
    Xm(:,t+1) = xm;
    Em(:,t+1) = xe - xm;
    Am(:,t) = am;

    %
    %Residual Norm Steepest Descend
    %
    vr = A'*rr;
    wr = A*vr;  ar = (vr'*vr)/(wr'*wr);
    xr = xr + ar*vr;
    rr = rr - ar*wr;

    %updating
    Rr(:,t+1) = rr;
    Vr(:,t+1) = vr;
    Wr(:,t+1) = wr;
    Xr(:,t+1) = xr;
    Er(:,t+1) = xe - xr;
    Ar(:,t) = ar;

end

% figure(1)
% title('Norm of the error')
% semilogy(vecnorm(Es))
% hold all
% semilogy(vecnorm(Em))
% semilogy(vecnorm(Er))
% legend('SD','MRI','RNSD')

%
% Minimized Functions
%
Xmin = min(min(Xs,Xm),Xr);
Xmax = max(max(Xs,Xm),Xr);

xmin = min(Xmin(1,:))-0.2;
xmax = max(Xmax(1,:))+0.2;
if (xmax - xmin)>(4*norm(xe))
    xmin = -3*norm(xe);
    xmax = 3*norm(xe);
end
ymin = min(Xmin(2,:))-0.2;
ymax = max(Xmax(2,:))+0.2;
if (ymax - ymin)>(4*norm(xe))
    ymin = -3*norm(xe);
    ymax = 3*norm(xe);
end
[X,Y] = meshgrid(xmin-0.1:0.1:xmax+0.1,ymin-0.1:0.1:ymax+0.1);
siz = size(X);
Fs = zeros(siz);
Fm = zeros(siz);

%F which minimizes Steepest Descend
for i=1:siz(1)
    for j=1:siz(2)
        Fs(i,j) = (A*([X(i,j);Y(i,j)]-xe))'*([X(i,j);Y(i,j)]-xe);
        Fm(i,j) = (b - A*([X(i,j);Y(i,j)]))'*(b - A*([X(i,j);Y(i,j)]));
    end
end

%
% Plot
%
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [0 0 20 25]);
figure(2)
% plot of X0
plot([ 0 Xs(1,1)], [0 Xs(2,1)],'b');
hold on
figure(3)
plot([ 0 Xm(1,1)], [0 Xm(2,1)],'r');
hold on
% figure(4)
plot([ 0 Xr(1,1)], [0 Xr(2,1)],'g');
hold on

%plot of Minimized Functions
figure(2)
cFs = [0.7 0.7 0.7];
contour(X,Y,Fs,'--','Levels',max(max(Fs))/50)
colormap(cFs)
axis equal
figure(3)
cFm = [0.3 0.3 0.3];
contour(X,Y,Fm,'-.','Levels',max(max(Fm))/50)
colormap(cFm)
axis equal
% figure(4)
% cFm = [0.3 0.3 0.3];
% contour(X,Y,Fm,'-.')
% colormap(cFm)
% axis equal


% annotation('textbox',[.9 .5 .4 .2], ...
%     'String','gray -- , F of SD, gray -., F of MRI and RNSD','EdgeColor','none')

figure(2)
grid on
axis([xmin-0.2 xmax+0.2 ymin-0.2 ymax+0.2])
title('Minimized Function - Steepest Descent')

figure(3)
legend({'MRI','RNSD'},'AutoUpdate','off','Location','northwest')
grid on
axis([xmin-0.2 xmax+0.2 ymin-0.2 ymax+0.2])
title('Min F - Minimal Residual Iteration - Residual Norm SD')

% figure(4)
% % legend({'SD'},'AutoUpdate','off','Location','northwest')
% grid on
% Xmin = min(min(Xs,Xm),Xr);
% Xmax = max(max(Xs,Xm),Xr);
% axis([min(Xmin(1,:))-0.2 max(Xmax(1,:))+0.2 min(Xmin(2,:))-0.2 max(Xmax(2,:))+0.2])
% title('Residual Norm Steepest Descend')


ARs = As.*Rs;
ARm = Am.*Rm;
ARr = Ar.*Rr;

for t=1:ni
    
    %Steepest Descend
    figure(2)
    plot([ Xs(1,t) Xs(1,t+1)], [Xs(2,t) Xs(2,t+1)],'b');
    
    %Minimal Residual Iteration
    figure(3)
    plot([ Xm(1,t) Xm(1,t+1)], [Xm(2,t) Xm(2,t+1)],'r');

%     figure(4)
    %Residual Norm Steepest Descend
    plot([ Xr(1,t) Xr(1,t+1)], [Xr(2,t) Xr(2,t+1)],'g');

end





















