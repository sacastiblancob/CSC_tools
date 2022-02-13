% %
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAM-SCHMIDT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Arnoldi-Modified Gram-Schmidt
% h = [1 3 1;1 3 1;2 6 2;2 6 2;3 8 3;3 8 3];
% A = diag(h(1:5,1),-1) + diag(h(:,2),0) + diag(h(1:5,3),1);
% v = [1;2;3;4;5;6];
% v = v./norm(v);
% m = 3;
% h = zeros(m,m);
% Q = zeros(length(A),m);
% Q(:,1) = v;
% for j=1:m
%     w = A*Q(:,j);
%     w0 = w;     %for comparing when reorthogonalization
%     for i=1:j
%         h(i,j) = w'*Q(:,i);
%         w = w - h(i,j)*Q(:,i);
%     end
%     if abs(norm(w) - norm(w0)) < eps(1E6) %reorthogonalization if needed (almost never needed)
%         w = -w;
%         for i=1:j
%             h(i,j) = w'*Q(:,i);
%             w = w - h(i,j)*Q(:,i);
%         end
%     end
%     if j~=m
%         h(j+1,j) = norm(w);
%         if h(j+1,j) <= 2*eps(1)
%             break
%         end
%    
%         Q(:,j+1) = w/h(j+1,j);
%     end
% end
% 
% %Proof (Figure 6.1 , page 181 Yousef Saad)
% AQ = A*Q;
% QH = Q*h;
% QH(:,m) = QH(:,m) + w;
% 
% % K = zeros(size(A));
% % K(:,1) = v/norm(v);
% % for i=2:n
% %     K(:,i) = A*K(:,i-1);
% %     K(:,i) = K(:,i)/norm(K(:,i));
% % end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HOUSEHOLDER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Mirror vector for Householder
% x = [1 2 3 4 5 6];
% x = x';
% [w,beta] = householderv(x);
% P = eye(length(x)) - beta*(w*w');
% z = P*x;

% %Householder Arnoldi (Yousef Saad)
% h = [1 3 1;1 3 1;2 6 2;2 6 2;3 8 3;3 8 3];
% A = diag(h(1:5,1),-1) + diag(h(:,2),0) + diag(h(1:5,3),1);
% v = [1;2;3;4;5;6];
% v = v/norm(v);
% z = v;
% n = length(v);
% m = 6;
% h = zeros(n,m+1);
% Q = zeros(n,m+1);
% 
% %naive
% I = eye(n);
% Pd = I;
% Pi = I;
% for j=1:m+1
%     w = zeros(n,1);
%     [w(j:end),beta] = householderv(z(j:end));
%     P = I - beta*(w*w');
%     h(:,j) = P*z;
%     Pd = Pd*P;
%     Pi = P*Pi;
%     Q(:,j) = Pd(:,j);
%     if j<=m
%         z = Pi*(A*Q(:,j));
%     end
% end
% 
% H = h(:,2:end);
% 
% K = zeros(size(A));
% K(:,1) = v/norm(v);
% for i=2:n
%     K(:,i) = A*K(:,i-1);
%     K(:,i) = K(:,i)/norm(K(:,i));
% end
% 
% %Proof (Figure 6.1 , page 181 Yousef Saad)
% AQ = A*Q;
% QH = Q(:,1:6)*H;

% %Arnoldi-Householder (Walker 1988) %NAIVE
% h = [1 3 1;1 3 1;2 6 2;2 6 2;3 8 3;3 8 3];
% A = diag(h(1:5,1),-1) + diag(h(:,2),0) + diag(h(1:5,3),1);
% v = [1;2;3;4;5;6];
% v = v/norm(v);
% z = v;
% n = length(v);
% m = 6;
% H = zeros(n,m+1);
% Q = zeros(n,m);
% I = eye(n);
% 
% [w,beta] = householderv(z);
% P = I - beta*(w*w');
% Pd = I;
% Pi = I;
% for j=1:m+1
%     Pd = Pd*P;
%     Pi = P*Pi;
%     H(:,j) = P*z;
% %     P
%     if j<=m
%         Q(:,j) = Pd(:,j);
%         z = Pi*(A*Q(:,j));
%         w = zeros(n,1);
%         [w(j+1:end),beta] = householderv(z(j+1:end));
%         P = I - beta*(w*w');
%     end
% end
% 
% %Proof for the case of Householder_Arnoldi
% % H = [h0,h1,h2,...,hm]   -> n x (m+1)
% % Q = [v1,v2,v3,...,vm]   -> n x m
% %
% % Hessemberg matrix of Arnoldi's
% % He = [h1,h2,h3,...,hm]  -> n x m
% % Proof is: A*Q = Q*He
% % This proof is only visible when m = n !!!
% % Because by the way it is storing things, for m < n a little modification
% % must be done, in orther to store all needed information
% %
% He = H(1:m,2:end);
% AQ = A*Q;
% QH = Q*He;
% A2 = Q*He*Q';

% % %Arnoldi-Householder (Walker 1988) %LITE
% h = [1 3 1;1 3 1;2 6 2;2 6 2;3 8 3;3 8 3];
% A = diag(h(1:5,1),-1) + diag(h(:,2),0) + diag(h(1:5,3),1);
% v = [1;2;3;4;5;6];
% % v = v/norm(v);
% z = v;
% n = length(v);
% m = 3;
% H = zeros(m+1,m+1);
% Q = zeros(n,m);
% 
% W = zeros(n,m+1);
% b = zeros(1,m+1);
% [W(:,1),b(1)] = householderv(z);
% for j=1:m+1
% 
%     %Computing H
%     H(1:j-1,j) = z(1:j-1);
%     
%     if j<=m  %this if presumibly can be removed, if Q dimensions are (n,m+1)
%         
%         H(j,j) = z(j) - b(j)*(z(j:n)'*W(j:n,j));
%         %Computing v_j, stored in j_th column of Q
%         % vj = P1*P2*...*PJ*ej
%         % i.e. succesive outer products by right side
%         q = zeros(n,1);
%         q(j) = 1;
%         for k = j:-1:1
%             bqTw = b(k)*(q(k:n)'*W(k:n,k));
%             for i = k:n
%                 Q(i,j) = q(i) - W(i,k)*bqTw;
%             end
%             q = Q(:,j);
%         end
% 
%         %Computing z_j+1
%         % z_j+1 = Pj*Pj-1*...*P2*P1*(A*v_j)
%         % i.e succesive outer products by left side
%         zaux = A*Q(:,j);
%         for k = 1:j
%             bzTw = b(k)*(zaux(k:n)'*W(k:n,k));
%             for i = k:n
%                 z(i) = zaux(i) - W(i,k)*bzTw;
%             end
%             zaux = z;
%         end
% 
%         %Computing new Householder vector (W(j+1)) and new beta
%         [W(j+1:end,j+1),b(j+1)] = householderv(z(j+1:end));
% 
%     end
% end
% 
% %Hessemberg Matrix (He) similar to A through Q*He*Q' when computed
% %completely, i.e., when m=n
% He = H(1:m,2:m+1);
% 
% %Proof for the case of Householder_Arnoldi
% % H = [h0,h1,h2,...,hm]   -> n x (m+1)
% % Q = [v1,v2,v3,...,vm]   -> n x m
% %
% % Hessemberg matrix of Arnoldi's
% % He = [h1,h2,h3,...,hm]  -> n x m
% % Proof is: A*Q = Q*He
% % This proof is only visible when m = n !!!
% % Because by the way it is storing things, for m < n a little modification
% % must be done, in orther to store all needed information
% %
% AQ = A*Q;
% QH = Q*He;
% A2 = Q*He*Q';
% 
% %
% % Outer products and related things
% %
% v = [1;2;3;4;5;6];
% v = v/norm(v);
% z = v;
% t = ones(length(z),1);
% [w,beta] = householderv(z);
% P = eye(length(z)) - beta*(w*w');
% Pt = zeros(length(z),1);
% btTw = beta*(t'*w);
% for i=1:length(z)
%     Pt(i) = t(i) - w(i)*btTw;
% end
% PT = P*t;
    
% Arnoldi Householder with optimized storage
h = [1 3 1;1 3 1;2 6 2;2 6 2;3 8 3;3 8 3];
A = diag(h(1:5,1),-1) + diag(h(:,2),0) + diag(h(1:5,3),1);
[Av,Ar,Ac] = full2csc(A);
v = [1;2;3;4;5;6];
m = 3;
[WH,b] = csc_arnoldi_householder(Av,Ar,Ac,v,m);
% n = length(v);
% z = v;
% WH = zeros(n+1,m+1);
% b = zeros(1,m+1);
% [WH(2:n+1,1),b(1)] = householderv(z);
% 
% for j=1:m+1
% 
%     %Computing H
%     WH(1:j-1,j) = z(1:j-1);
%     WH(j,j) = z(j) - b(j)*(z(j:n)'*WH(j+1:n+1,j));
%     
%     if j<=m  %this if presumibly can be removed, if Q dimensions are (n,m+1)
%         
%         %Computing v_j, stored in j_th column of Q
%         % vj = P1*P2*...*PJ*ej
%         % i.e. succesive outer products by right side
%         v = zeros(n,1);
%         q = v;
%         q(j) = 1;
%         for k = j:-1:1
%             bqTw = b(k)*(q(k:n)'*WH(k+1:n+1,k));
%             for i = k:n
%                 v(i) = q(i) - WH(i+1,k)*bqTw;
%             end
%             q = v;
%         end
% 
%         %Computing z_j+1
%         % z_j+1 = Pj*Pj-1*...*P2*P1*(A*v_j)
%         % i.e succesive outer products by left side
%         zaux = csc_matvec(Av,Ar,Ac,v);
%         for k = 1:j
%             bzTw = b(k)*(zaux(k:n)'*WH(k+1:n+1,k));
%             for i = k:n
%                 z(i) = zaux(i) - WH(i+1,k)*bzTw;
%             end
%             zaux = z;
%         end
% 
%         %Computing new Householder vector (W(j+1)) and new beta
%         [WH(j+2:n+1,j+1),b(j+1)] = householderv(z(j+1:end));
%     end
% end

















