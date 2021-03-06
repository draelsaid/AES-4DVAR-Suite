function [kappaSx,kappaSp,kappaD] = HessNumCondNoCalc(n,N,fcD,fcR,fcM,obs_f_x,mlvl)

% Cost function matrices
O = zeros(N,N); I = eye(N);

% D Matrix
[Dinv{1:n+1,1:n+1}] = deal(O); [rtD{1:n+1,1:n+1}] = deal(O);

[~,~,B,Q] = fcD(randn(N,1));
[~,~,Rinv] = fcR(rand(N,n+1)); sigma_o = (1/Rinv(1,1))^0.5;

Binv = inv(B);
Qinv = inv(Q);

Dinv{1,1} = Binv;
rtD{1,1} = B^0.5;

for i=2:n+1
      Dinv{i,i} = Qinv;
       rtD{i,i} = Q^0.5;
end
rtD = cell2mat(rtD);
Dinv = cell2mat(Dinv);

% %%%% FOR ADVECTION EQUATION %%%%%%

M = fcM()^mlvl;            % Model in matrix form (must be LINEAR OR LINEARISED)

%%%%%%%%%%%%%%%%%%%%%%%%% L, Linv, H, H^T %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[L{1:n+1,1:n+1}] = deal(O);
for i=1:n+1
    L{i,i} = I;
end

for i=2:n+1
L{i,i-1} = -M;
end

L = cell2mat(L);

[Linv{1:n+1,1:n+1}] = deal(O);
for i=1:n+1
  for j=1:i
    Linv{i,j} = M^(i-j);
  end
end
Linv = cell2mat(Linv);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % %%%%%% FOR Model matrix of L95 (non-linear) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fnl = @l95; fnlmtx = @l95mtx; ftl = @l95tlmtx; schtl = @rk4tlmtx; schnl = @rk4nl; schnlmtx = @rk4nlmtx; 
% [M{1,1:n}] = deal(O); xl(:,1) = xi; ddx(:,1) = rand(N,1)*0.01;

%NON LIN%
% for i=1:n 
% xl(:,i+1) = Mnl(xl(:,i),fnlmtx,schnlmtx,dt,1);    
% [xl(:,i+1),~,Mm] = Mnlmtx(xl(:,i),fnlmtx,schnlmtx,dt,1);
%               M{1,i} = Mm;
% end

%LINEARISED%
% for i=1:n 
% xl(:,i+1) = Mnl(xl(:,i),fnl,schnl,dt,1);    
% [ddx(:,i+1),~,Mm] = Mtlmtx(ddx(:,i),xl(:,i),ftl,fnl,schtl,schnl,dt,1);
%               M{1,i} = Mm;
% end

% %%%%%%%%%%%%%%%%%%%%%%%% L, Linv, H, H^T %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [L{1:n+1,1:n+1}] = deal(O);
% for i=1:n+1
%     L{i,i} = I;
% end
% 
% for i=2:n+1
% L{i,i-1} = -M{1,i-1};
% end
% 
% L = cell2mat(L);
% 
% [Linv{1:n+1,1:n+1}] = deal(O);
% 
% for i=1:n+1
%    Linv{i,i} = I;
% end
% 
% for k=1:n
%  for i=1:n-k+1
%      Mdum=I;
%    for j=i:k+i-1
%           Mdum = M{1,j}*Mdum;
%    end      
%     Linv{i+k,i} = Mdum;
%   end
% end
% 
% for i=1:n+1
%     condMM(i) = cond(Linv{i,1});
% end
% 
% Linv = cell2mat(Linv);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=zeros(obs_f_x,N); j=1; 
  for i=1:obs_f_x:N
    h(j,i) = 1;
    j=j+1;
  end

[H{1:n+1,1:n+1}] = deal(zeros(size(h)));
  for j=1:n+1
    H{j,j} = h;
  end
H = cell2mat(H);

%% Hessians
G = (Linv'*H')*H*Linv;

Sp = Dinv + sigma_o^(-2)*G;                           % p formulation Hessian
%PSp = eye(size(H'*H)) + sigma_o^(-2)*(rtD*G*rtD);    % precond p formulation Hessian
Sx = L'*Dinv*L + sigma_o^(-2)*(H'*H);                 % x formulation Hessian

% Condition number calculations
kappaSp = cond(Sp);
%kappaPSp = cond(PSp);
kappaSx = cond(Sx);
kappaD = cond(Dinv);
%kappaLTDL = cond(L'*Dinv*L);
%kappaSpterm = cond(G);
end

