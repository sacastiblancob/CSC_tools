%SOLVER PARA UN DOMINIO COMPUTACIONAL 2D DE LA ECUACION DIFERENCIAL DE
%TRANSPORTE - Metodo pasos fraccionados (PF)
%
%DOMINIO ESPACIAL   -a < x < b
%                   -a < y < b
%
%DOMINIO TEMPORAL   tincial < t < tfinal
%   (donde t inicial es ligeramente mayor a t0)
%
%   En este caso particular se toma t0 = 1s (ver funcion analitica), y
%   tinicial = 1s + dt
%
%Matriz de rigidez 'optimizada', sacando de ella las condiciones de
%frontera, e instertando estas en el vector de valores fijos Caux

Dx=0.25;
Dy=0.25;
u=0;
v=0;

%Discretizacion espacial
a=4;
b=4;
dx=0.25;
dy=0.25;

n=round((a-(-b))/dx+1);
m=round((a-(-b))/dy+1);
x=-a:dx:b;
y=-a:dy:b;
X=kron(ones(1,m),x);
Y1=kron(y,ones(1,n));
Y=Y1(n*m:-1:1);
MX=reshape(X,m,n)';
MY=reshape(Y,m,n)';

%Discretizacion temporal
dt=.025;
tinicial=1+dt;
tfinal=10;
t=tinicial+dt:dt:tfinal;

%Constantes de la ecuacion diferencial computacional
Sx=(Dx*dt)/(dx^2);
Sy=(Dy*dt)/(dy^2);
CFLx=(u*dt)/dx;
CFLy=(v*dt)/dy;
CFL = max(CFLx,CFLy);
%Comprobacion de los parametros Sx, Sy, CFLx y CFLy
%     if CFL <=1
%         disp("Los valores de Sx, Sy, CFLx y CFLy estan dentro de los limites adecuados")
%     else
%         disp("ATENCION: Los valores de Sx, Sy, CFLx y CFLy no estan dentro de los limites adecuados")
%         return
%     end

%Valores de los nodos que estan en la frontera
front=[1:n,n+1:n:n*m-2*n+1,2*n:n:n*m-n,n*m-n+1:n*m];
front=sort(front);
fronts=(n+3:2*n-2);
frontd=(3*n-1:n:n*m-2*n-1);
frontin=(n*m-2*n+3:n*m-n-2);
frontiz=(2*n+2:n:n*m-3*n+2);

%Ensamblaje matriz sistema de ecuaciones para el paso 2 del planteamiento
%numerico
ap=(1+2*Sx+2*Sy);
ax=-Sx;
ay=-Sy;

% % % K full     
% K11=diag(ones(n-2,1));
% K12=ax*diag(ones(n-3,1),1)+ax*diag(ones(n-3,1),-1)+ap*diag(ones(n-2,1));
% K21=diag(ones(n-3,1),1)+diag(ones(n-3,1),-1);
% K22=ay*diag(ones(n-2,1));
% K=kron(K11,K12)+kron(K21,K22);
% % % dimen(iv) = length(K);
% %CondK=cond(K);
% %disp('Condicion de la matriz K es: '); disp(CondK);
% K=sparse(K);
% diagk = sparse(diag(diag(K)));

% K CSC
    [K11v,K11r,K11c] = csc_diag(ones(n-2,1),0);
    h1 = [ax*ones(n-2,1) ap*ones(n-2,1) ax*ones(n-2,1)];
    d1 = [-1 0 1];
    [K12v,K12r,K12c] = csc_diag(h1,d1);
    h2 = [ones(n-3,1) ones(n-3,1)];
    d2 = [-1 1];
    [K21v,K21r,K21c] = csc_diag(h2,d2);
    [K22v,K22r,K22c] = csc_diag(ay*ones(n-2,1),0);

    [K1v,K1r,K1c] = csc_kron(K11v,K11r,K11c,K12v,K12r,K12c);
    [K2v,K2r,K2c] = csc_kron(K21v,K21r,K21c,K22v,K22r,K22c);

    [Kv,Kr,Kc] = csc_sum(K1v,K1r,K1c,K2v,K2r,K2c);
    clearvars K11v K11r K11c K12v K12r K12c K21v K21r K21c K22v K22r K22c K1v K1r K1c K2v K2r K2c
%     dimen(iv) = length(Kc)-1;

% tic

% %Jacobi
[Pv,Qv,Qc,Qr] = csc_prejacobi(Kv,Kr,Kc);

% %SOR
% [Pv,Pr,Pc,Qv,Qr,Qc] = csc_preSOR(Kv,Kc,Kr,1.1);

%     %Computing packed LU decomposition for matrix K
%     LU = packlu(K);
% %     [Kv,Kr,Kc] = full2csc(K);
% %     [LUv,LUr,LUc] = packlu_csc(Kv,Kr,Kc);
% tocsit(iv) = toc;
%     toc

%spy(K)

%Solver

Co=analitica(X,Y,tinicial,u,-v,Dx,Dy);
conte=1;
sizet=size(t);
errorPF=zeros(1,sizet(2));
Caux=ones(1,m*n-2*m-2*n+4);
trash=zeros(n*m);
times = 0;

for i=tinicial+dt:dt:tfinal
    Ca=analitica(X,Y,i,u,-v,Dx,Dy);
    l=1;
    for j=1:m*n
        if find(j==front)
            trash(j)=0;
        elseif find(j==fronts)
            Caux(l)=Co(j-1)*CFLx + Co(j-n)*CFLy + Co(j)*(1-CFLx-CFLy) - ay*Ca(j-n);
            l=l+1;
        elseif find(j==frontd)
            Caux(l)=Co(j-1)*CFLx + Co(j-n)*CFLy + Co(j)*(1-CFLx-CFLy) - ax*Ca(j+1);
            l=l+1;
        elseif find(j==frontin)
         Caux(l)=Co(j-1)*CFLx + Co(j-n)*CFLy + Co(j)*(1-CFLx-CFLy) - ay*Ca(j+n);
            l=l+1;
        elseif find(j==frontiz)
            Caux(l)=Co(j-1)*CFLx + Co(j-n)*CFLy + Co(j)*(1-CFLx-CFLy) - ax*Ca(j-1);
            l=l+1;
        elseif j==n+2
            Caux(l)=Co(j-1)*CFLx + Co(j-n)*CFLy + Co(j)*(1-CFLx-CFLy) - ay*Ca(j-n)-ax*Ca(j-1);
            l=l+1;
        elseif j==2*n-1
            Caux(l)=Co(j-1)*CFLx + Co(j-n)*CFLy + Co(j)*(1-CFLx-CFLy) - ay*Ca(j-n)-ax*Ca(j+1);
            l=l+1;
        elseif j==n*m-2*n+2
            Caux(l)=Co(j-1)*CFLx + Co(j-n)*CFLy + Co(j)*(1-CFLx-CFLy) - ay*Ca(j+n)-ax*Ca(j-1);
            l=l+1;
        elseif j==n*m-n-1
            Caux(l)=Co(j-1)*CFLx + Co(j-n)*CFLy + Co(j)*(1-CFLx-CFLy) - ay*Ca(j+n)-ax*Ca(j+1);
            l=l+1;
        else
            Caux(l)=Co(j-1)*CFLx + Co(j-n)*CFLy + Co(j)*(1-CFLx-CFLy);
            l=l+1;
        end
    end

    if i == tinicial+dt
        C = Caux';
    end

    tol = 1E-6;
    tic

    % Matlab method
    % C = K\(Caux');

    % Matlab pcg
    % [x,flag,relres,t] = pcg(K,Caux',tol,100,diagk,speye(length(K)),C);

    % Gaussian elimination
    % C = elgauss(K,Caux');

    % Pack-LU
    % C = solpacklu(LU,Caux');

    % Pack-LU CSC
    % C = solpacklu_csc(LUv,LUr,LUc,Caux');

    % %Jacobi CSC
    [C,t] = csc_jacobi(Kv,Kr,Kc,Pv,Qv,Qc,Qr,Caux',C,100,tol);

    % %SOR CSC
    % [C,t] = csc_SOR(Kv,Kr,Kc,Pv,Pr,Pc,Qv,Qr,Qc,Caux',C,100,tol);

    % % CG CSC
    % [C,t] = csc_CG(Kv,Kr,Kc,Caux',C,100,tol);

    % % DPCG CSC
    % [C,t] = csc_DPCG(Kv,Kr,Kc,Caux',C,100,tol);

    times = times + 1;

    CM = reshape(C,n-2,m-2);
    C1 = [pi*ones(1,n-2);CM;pi*ones(1,n-2)];
    C1 = C1(:)';
    C2 = [pi*ones(1,n),C1,pi*ones(1,n)];
    C2(C2==pi)=Ca(C2==pi);

%         %MCc=reshape(C,m,n)';
%         MCc = reshape(C2,m,n)';
% %         MCa=reshape(Co,m,n)'; 
%         surf(MX,MY,MCc); axis([-4 4 -4 4 -1 8]); view(2)
%         drawnow;
%         %pause(dt)
%         %errorPF(conte)=norm(C2-Ca);
%         %conte=conte+1;
    Co=C2;
end

%toc
%ultimos valores de concentracion calculados reorganizados en matrices

%     MCc=reshape(C2,m,n)';                            %C calculada reordenada en matriz 
%     MCa=reshape(Ca,m,n)';                           %C analitica reordenada en matriz
%     MdifPF=reshape(Ca-C2,m,n)';                        %Diferencia entra C analitica y C calculada reordenada en matriz
% 
%     % plot(t,log10(errorPF))                          %Grafica de la evolucion del error en el tiempo
%     % xlabel('Tiempo (t)')
%     % ylabel('Log10(error)')
% 
%     conv=log10(dt+dx^2);                             %Orden de convergencia
%     faccon=[Sx,Sy,CFLx,CFLy,Sx+Sy,CFLx+CFLy];        %Factores de convergencia

%end