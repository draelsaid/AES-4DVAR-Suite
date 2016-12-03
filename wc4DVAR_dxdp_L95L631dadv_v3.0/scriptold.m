% Weak Constraint 4DVAR main script
% Adam El-Said, University of Reading (06/2013)

clear all
% Close any open figure windows
g=gcf;
close all

%% Graph position settings
ss=get(0,'ScreenSize'); fig_width=0.55*ss(3); fig_height=0.4*ss(4);
pos_1=[0.02*ss(3),0.5*ss(4),fig_width,fig_height]; pos_2=[0.02*ss(3),0.05*ss(4),fig_width,fig_height]; pos_3=[0.6*ss(3),0.45*ss(4),0.4*ss(3),0.45*ss(4)];
pl=0.6*ss(3); pbb=0.05*ss(4); pw=0.3*ss(3); ph=0.35*ss(4);
pos_4=[pl,pbb,pw,ph]; pos_5=[0.02*pw,0.8*ph,0.95*pw,0.2*ph]; pos_6=[0.02*pw,0.02*ph,0.95*pw,0.75*ph];

% Model settings
N=100; n=60; % Size of spatial and temporal domains
dt = 0.01;    % 0.025 dt => 3 hour time step (L95), 1 time unit is 5 days
dx = 1/N;     % Spatial interval size
C = (-1*dt)/(dx)

%% TRUTH
fnlt = @adv1d; scht = @lpde;    % TRUE model and scheme func handles

% IF THE MODEL IS LINEAR, CHECK MINIMISATION PROCEDURE SUITABILITY.

% True Covariance Matrices D and R
sigmabt = 0.002; Lsbt = 0.03;    % True sigma_b and lengthscale for C_B
sigmaqt = 0.001; Lsqt = 0.02;    % True sigma_q and lengthscale for C_Q
sigmaot = 0.001;                 % True sigma_o
fcBt = @(N) soar(N,Lsbt,dx);                  % B MATRIX function handle
fcQt = @(N) laplacian(N,Lsqt,dx);             % Q MATRIX function handle
fDt = @(x) dmtx(x,sigmabt,sigmaqt,fcBt,fcQt); % D MATRIX func handle
fRt = @(x) rmtx(x,sigmaot);                   % R MATRIX function handle

%% Initial Conditions

% Gaussian-shaped curve (for 1D advection equation ONLY)
xc = (N*dx)/2; sigma = 0.2; Ampl = 2; % Mode, Standard dev and Amplitude
xs = linspace(0,1,N);                 % Physical space index 

for i=1:N
    xto(i) = Ampl*exp(-(xs(i)-xc).^2/sigma^2);
end

% % For chaotic ODE models
% xto = Mnl(rand(N,1),fnlt,scht,dt,50); % Model 'spin up' for initial state

% Truth generation
[~,etat] =  fDt(rand(N,n+1));       % True model errors
bias = 0; etat = bias + etat;   % Model error bias
etat(:,1) = 0;                      % Model error time adj (eta_0 not exist)

% Truth model trajectory with model error
xttraj(:,1) = xto;
for i=1:n
xttraj(:,i+1) = Mnl(xttraj(:,i),fnlt,scht,dt,n) + etat(:,i+1);
end

etat = etat(:,2:end); % Model error time dimension adjustment

% % Time-movie for advection equation
% for i=1:n
%  figure(10), plot(xs,xttraj(:,i),'b'), title('Truth run'), xlabel('Position'), legend('Analysis from t_0 to t_{now}'), ylabel('Value')
%  drawnow
% end

%% ASSIMILATION

% Model function handle and integration scheme setup
fnl = @adv1d;      schnl = @lpde;
ftl = @adv1d;      schtl = @lpde;
fadj= @adv1d_adj;  schadj = @lpde_adj;

% Model propogators
fnlm = @(x) Mnl(x,fnl,schnl,dt,1);
ftlm = @(dx,x) Mtl(dx,x,ftl,fnl,schtl,schnl,dt,1); 
fadjm = @(dy,y) Madj(dy,y,fadj,fnl,schadj,schnl,dt,1);

fnlmsc = @(x) Mnl(x,fnl,schnl,dt,n);

% Background error creation [truth plus some noise ~N(0,D)]
[~,xbe] = fDt(rand(N,n+1));
pb = xttraj;% + xbe;
pb(:,2:end)=0;             % Background model errors eta^b=0

[~,xbtraj] = Mnl(pb(:,1),fnlt,scht,dt,n);

%% Observation creation parameters

[~,ye] = fRt(rand(N,n+1));    
y = xttraj;% + ye;           % Obs = truth + noise ~N(0,R)

obs_f_t = 1;      % Obs frequency (time) (Regular spacing. =1 obs all times)
obs_f_x = 1;      % Obs frequency (space) (Regular spacing. =1 full obs)
obsxo = 1;        % Determines observation starting point (space)
obsto = 1;        % Determines observation starting point (time)

fx = floor(N/obs_f_x); % Variable sizes determined by input

fH = @(y) Hop(y,obs_f_x,obs_f_t,obsxo,obsto);         % H operator
fHT = @(y) HopT(y,N,n+1,obs_f_x,obs_f_t,obsxo,obsto); % H operator Transposed

y = fH(y);

%% Covariance Matrices (ASSIMILATION)

sigmab = 0.2;  Lsb = 0.03;    % Estimated sigmab and lengthscale for B
sigmaq = 0.001;  Lsq = 0.02;  % Estimated sigmaq and lengthscale for Q
sigmao = 0.001;               % Estimated sigmao

fcB = @(N) soar(N,Lsb,dx);         % C_B MATRIX function handle
fcQ = @(N) laplacian(N,Lsq,dx);    % C_Q MATRIX function handle

fD = @(x) dmtx(x,sigmab,sigmaq,fcB,fcQ); % D MATRIX function handle
fR = @(x) rmtx(x,sigmao);                % R MATRIX function handle

%% Gradient Test

% % Function handles
%       cfdp = @(dp) wcdp(dp,pb,y,fH,fHT,fD,fR,fnlm,fadjm);
%       cfdx = @(dx) wcdx(dx,pb,y,fH,fHT,fD,fR,fnlm,fadjm);
%  cfsc4dvar = @(dw) sc4dvar(dw,pb,y,fH,fHT,fD,fR,fnlmsc,fadjm);
   
% % Tests
%        phidp = testgrad(cfdp,pb); % Tests wc4dvar dp formulation
%        phidx = testgrad(cfdx,pb); % Tests wc4dvar dx formulation
%  phisc4dvar = testgrad(cfsc4dvar,pb,1); % Tests sc4dvar

%% Minimisation

% Minimisation parameters
xstr='Maximum number of iterations'; ystr='Tolerance'; zstr={'100','1d-4'};
min_str=inputdlg({xstr,ystr},'Minimization control',1,zstr); 
max_iterations=str2num(min_str{1}); tolerance=str2num(min_str{2});

p = pb; x = pb;            % wc4DVAR Initial guess
xo = p(:,1);               % sc4DVAR initial guess

s = size(pb,2);
p = reshape(p,s*N,1);
x = reshape(x,s*N,1);

%% LINEAR CONJUGATE GRADIENT (requires variable order [DJ,J] in 4dvar)

% Minimsation for sc4DVAR
[bbsc,~] = sc4dvar(zeros(size(xo)),pb,y,fH,fHT,fD,fR,fnlmsc,fadjm); bbsc=-bbsc;
scfun = @(xo) sc4dvar(xo,pb,y,fH,fHT,fD,fR,fnlmsc,fadjm) + bbsc;

[X_asc,FLAGsc,relressc,~,dJXsc,xmatsc] = minimise_pcg(scfun,bbsc,tolerance,max_iterations,[],[],xo);

% sc4DVAR Solution
xatrajsc = Mnlsc(X_asc,fnlm,n);

% Minimsation for dp
[bbdp,~] = wcdp(zeros(size(p)),pb,y,fH,fHT,fD,fR,fnlm,fadjm); bbdp=-bbdp;
wcdpfun = @(p) wcdp(p,pb,y,fH,fHT,fD,fR,fnlm,fadjm) + bbdp;

[X_adp,FLAGdp,relresdp,~,dJXdp,xmatdp] = minimise_pcg(wcdpfun,bbdp,tolerance,max_iterations,[],[],p);

% dp Solution
X_adp = reshape(X_adp,N,s);
xatrajdp = Lopnlinv(X_adp,fnlm);
eta = X_adp(:,2:end); 

% Minimisation for dx
[bbdx,~] = wcdx(zeros(size(x)),pb,y,fH,fHT,fD,fR,fnlm,fadjm); bbdx=-bbdx;
wcdxfun = @(x) wcdx(x,pb,y,fH,fHT,fD,fR,fnlm,fadjm) + bbdx;

[X_adx,FLAGdx,relresdx,~,dJXdx,xmatdx] = minimise_pcg(wcdxfun,bbdx,tolerance,max_iterations,[],[],x);

% dx Solution
 X_adx = reshape(X_adx,N,s);
 xatrajdx = X_adx;

% Fix iteration count and vector arrangement for plots
lll = size(dJXdp); itdp = lll(1);
ll = size(dJXdx); itdx = ll(1);
l = size(dJXsc); itsc = l(1);
dJXdp = dJXdp'; dJXdx = dJXdx'; dJXsc = dJXsc'; 

% %% MinRes minimisation

% % Minimsation for dp
% [bbdp,~] = wcdp(zeros(size(p)),pb,y,fH,fHT,fD,fR,fnlm,fadjm); bbdp=-bbdp;
%  wcdpfun = @(p) wcdp(p,pb,y,fH,fHT,fD,fR,fnlm,fadjm) + bbdp;
% 
% [X_adp,FLAGdp,relresdp,itdp,dJXdp,xmatdp] = minres(wcdpfun,bbdp,tolerance,max_iterations,[],[],p);

% % dp Solution
% X_adp = reshape(X_adp,N,s);
% xatrajdp = Lopnlinv(X_adp,fnlm);
% eta = X_adp(:,2:end); 

% % Minimsation for dx
% [bbdx,~] = wcdx(zeros(size(x)),pb,y,fH,fHT,fD,fR,fnlm,fadjm); bbdx=-bbdx;
%  wcdxfun = @(x) wcdx(x,pb,y,fH,fHT,fD,fR,fnlm,fadjm) + bbdx;
% 
% [X_adx,FLAGdx,relresdx,itdx,dJXdx,xmatdx] = minres(wcdxfun,bbdx,tolerance,max_iterations,[],[],x);

% % dx Solution
% X_adx = reshape(X_adx,N,s);
% xatrajdx = Lopnlinv(X_adx,fnlm);

%% Compute J (function value) for each iterate (convergence plot)
%  (This is needed for MinRes and min_pcg)
for i=1:size(xmatdp,2)
    qq=xmatdp(:,i);
    [~,Jdp]=wcdp(qq,pb,y,fH,fHT,fD,fR,fnlm,fadjm);
    JXdp(i)=Jdp;
end

for i=1:size(xmatdx,2)
    qq=xmatdx(:,i);
    [~,Jdx]=wcdx(qq,pb,y,fH,fHT,fD,fR,fnlm,fadjm);
    JXdx(i)=Jdx;
end

for i=1:size(xmatsc,2)
    qq=xmatsc(:,i);
    [~,Jsc]=sc4dvar(qq,pb,y,fH,fHT,fD,fR,fnlmsc,fadjm);
    JXsc(i)=Jsc;
end

%% POLAK REBIERRE NON-LINEAR CONJUGATE GRADIENT (requires [J,DJ] in 4dvar
% procedure IN THAT ORDER)

% % Minimisation for dp
% [X_adp JXdp dJXdp itdp] = minimize_mod_crit(p,'wcdp',max_iterations,tolerance,pb,y,fH,fHT,fD,fR,fnlm,fadjm);
% 
% X_adp = reshape(X_adp,N,s);
% xatrajdp = Lopnlinv(X_adp,fnlm);
% 
% % Minimisation for dx
% [X_adx JXdx dJXdx itdx] = minimize_mod_crit(x,'wcdx',max_iterations,tolerance,pb,y,fH,fHT,fD,fR,fnlm,fadjm);
% 
% X_adx = reshape(X_adx,N,s);
% xatrajdx = X_adx;

%% Trajectory time adjustment
% xttraj = xttraj(1:N,2:end);
% xbtraj = xbtraj(1:N,2:end);
% xatrajdx = xatrajdx(1:N,2:end);
% xatrajdp = xatrajdp(1:N,2:end);

%% Error calculations

% Relative total error
terrordx = norm(xttraj - xatrajdx) / norm(xttraj)
terrordp = norm(xttraj - xatrajdp) / norm(xttraj)
terrorsc = norm(xttraj - xatrajsc) / norm(xttraj)

% Relative error of each time interval
for i=1:s
errordx(i) = norm(xttraj(:,i) - xatrajdx(:,i))/norm(xttraj(:,i));
end

for i=1:s
errordp(i) = norm(xttraj(:,i) - xatrajdp(:,i))/norm(xttraj(:,i));
end

for i=1:s
errorsc(i) = norm(xttraj(:,i) - xatrajsc(:,i))/norm(xttraj(:,i));
end

% Relative absolute error spatially and temporally
for i=1:N
    for j=1:s
      cerrordx(i,j) = abs(xttraj(i,j) - xatrajdx(i,j))/abs(xttraj(i,j));
    end
end

for i=1:N
    for j=1:s
      cerrordp(i,j) = abs(xttraj(i,j) - xatrajdp(i,j))/abs(xttraj(i,j));
    end
end

for i=1:N
    for j=1:s
      cerrorsc(i,j) = abs(xttraj(i,j) - xatrajsc(i,j))/abs(xttraj(i,j));
    end
end

%% L95 plots

% % Contour plots of truth, background and solution
% contourf(xatrajdx'), title('Contour plot Lorenz95 (J(x) Solution)'), ylabel('Time (n)'), xlabel('Variables (X_i)'); 
% figure; contourf(xatrajdp'), title('Contour plot Lorenz95 (J(p) Solution)'), ylabel('Time (n)'), xlabel('Variables (X_i)');
% figure; contourf(xttraj'), title('Contour plot Lorenz95 (TRUTH)'), ylabel('Time (n)'), xlabel('Variables (X_i)');
% figure; contourf(xbtraj'), title('Contour plot Lorenz95 (Background)'), ylabel('Time (n)'), xlabel('Variables (X_i)'); 
% 
% % J(x), J(p) spatio-temporal contour error plots
% figure; 
% e1 = subplot(2,1,1);
% contourf(cerrordx), title('Contour plot Lorenz95 (J(x) errors)'), xlabel('Time (n)'), ylabel('Variables (X_i)'), colorbar;
% e2 = subplot(2,1,2);
% contourf(cerrordp), title('Contour plot Lorenz95 (J(p) errors)'), xlabel('Time (n)'), ylabel('Variables (X_i)'), colorbar;
% %set([e1 e2],'clim',[0 max( max(max(log(cerrordp))),  max(max(log(cerrordx))) )]);
% 
% % J(x), J(p) relative space-averaged contour error plots vs time 
% figure;
% plot(errordp), axis([-inf,inf,-inf,inf]), title('Relative error'), xlabel('Time (n)'), ylabel('Error value')
% hold; plot(errordx,'r'), plot(errorsc,'g'), legend('Errors J(p)', 'Errors J(x)', 'Errors J(x_0)'); hold;

%% Advection equation plots

% J(x_0), J(p), J(x) relative space-averaged contour error plots vs time 
figure;
plot(errordp,'b','LineWidth',1.5), axis([-inf,inf,-inf,inf]), title('Relative error'), xlabel('Time (n)'), ylabel('Error value')
hold; plot(errordx,':r','LineWidth',1.5), plot(errorsc,'g','LineWidth',1.5), ...
legend('Errors J(p)', 'Errors J(x)', 'Errors J(x_0)'); hold;

% xt, xa, xb time series
figure('Position',pos_1)
i = 0; ts = 0;
for k = 1:4
    if k==4 
        i=s; ts=i*dt;
    else
i=ceil((k-1)*(n/4) + 1); 
ts = (i-1)*dt;
    end
subplot(2,2,k), hold on, ...
title(['t=' num2str(ts)]), plot(xs,xatrajdp(:,i),'b','LineWidth',1.5), ...
plot(xs,xatrajdx(:,i),'--r','LineWidth',1.5), plot(xs,xatrajsc(:,i),'g','LineWidth',1.5), ...
plot(xs,xttraj(:,i),'k:','LineWidth',2.5), xlabel('x'), ylabel('Value'), ...
end
legend('Solution J(p)','Solution J(x)','Solution J(x_0)','Truth','Location',[0.86 0.85 0.1 0.1])

% eta^t and eta time series (model error comparison)
figure('Position',pos_3)
i = 0; ts = 0;
for k = 1:4
    if k==4 
i = n; ts = i*dt;
    else
i = ceil((k-1)*(n/4) + 1); 
ts = (i-1)*dt;
    end
subplot(2,2,k), hold on, plot(xs,etat(:,i)), plot(xs,eta(:,i),'r'), xlabel('x'), ...
ylabel('Value'), title(['\eta_{' num2str(i) '}']), ...
end
legend('\eta^t \times \Delta t','\eta','Location',[0.86 0.89 0.1 0.1])
% qbar = num2str(mean(mean((eta)))); qbart = num2str(mean(mean((etat)))); qbar = qbar(1:6); qbart = qbart(1:6);
% text(-0.5,-0.25,['b=' num2str(qbar) ', b^t=' num2str(qbart)],'sc')

%% Convergence plot J(p)
h3=figure('Position',pos_6); clf;
subplot(2,1,1)
semilogy(linspace(1,max(size(JXdp)),max(size(JXdp))),JXdp), title('Convergence of cost function, J(p)'), xlabel('Iteration'), ylabel('Cost function'), axis([1,inf,0,inf]),
subplot(2,1,2)
for i = 1 : itdp
  g =  dJXdp(:,i);
  graddp(i)=norm(g);
end
semilogy(linspace(1,max(size(graddp)),max(size(graddp))),graddp), title('Convergence of gradient, J(p)'), xlabel('Iteration'), ylabel('Norm of gradient'), axis([1,inf,0,inf]),

%% Convergence plot J(x)
h3=figure('Position',pos_6); clf; 
subplot(2,1,1)
semilogy(linspace(1,max(size(JXdx)),max(size(JXdx))),JXdx,'r'), title('Convergence of cost function, J(x)'), xlabel('Iteration'), ylabel('Cost function'), axis([1,inf,0,inf]),
subplot(2,1,2)
for i = 1 : itdx
  g =  dJXdx(:,i);
  graddx(i)=norm(g);
end
semilogy(linspace(1,max(size(graddx)),max(size(graddx))),graddx,'r'), title('Convergence of gradient, J(x)'), xlabel('Iteration'), ylabel('Norm of gradient'), axis([1,inf,0,inf]),

%% Convergence plot J(x_0)
h3=figure('Position',pos_6); clf; 
subplot(2,1,1)
semilogy(linspace(1,max(size(JXsc)),max(size(JXsc))),JXsc,'g'), title('Convergence of cost function, J(x_0)'), xlabel('Iteration'), ylabel('Cost function'), axis([1,inf,0,inf]),
subplot(2,1,2)
for i = 1 : itsc
  g =  dJXsc(:,i);
  gradsc(i)=norm(g);
end
semilogy(linspace(1,max(size(gradsc)),max(size(gradsc))),gradsc,'g'), title('Convergence of gradient, J(x_0)'), xlabel('Iteration'), ylabel('Norm of gradient'), axis([1,inf,0,inf]),
