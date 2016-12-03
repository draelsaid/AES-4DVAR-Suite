clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Weak-Constraint 4D-Var with model error forcing applied to 
%  the 1D Advection model
%
%  Written by Adam El-Said 2012
%
%% Assimilation variables
%    xb:         Background term
%    h:          Spatial component of observation operator
%    ht:         Temporal component of observation operator
%    B:          Background weighting matrix
%    R:          Observation weighting matrix
%    Q:          Model Error weighting matrix
%    [z]:        [x_0,eta_1,...,eta_N] analysis
%    [ystore]:   Observation values
%
%% Model variables
%    L:          Physical domain size
%    a:          Velocity parameter in advection equation
%    dt:         Time step size
%    dx:         Spatial step size
%    nt:         Number of time steps
%    nx:         Number of gridpoints
%
%% Output:
%    [J,DJ]: Cost function and gradient
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------
%  1. Graph positions, model parameters and vector allocations
%-------------------------------------------------------------------------
% rng('default') % Keeps the same noise used

% Close any open figure windows
g=gcf;
close all

%% Graph position settings
ss=get(0,'ScreenSize'); fig_width=0.55*ss(3); fig_height=0.4*ss(4);
pos_1=[0.02*ss(3),0.5*ss(4),fig_width,fig_height]; pos_2=[0.02*ss(3),0.05*ss(4),fig_width,fig_height]; pos_3=[0.6*ss(3),0.45*ss(4),0.4*ss(3),0.45*ss(4)];
pl=0.6*ss(3); pb=0.05*ss(4); pw=0.3*ss(3); ph=0.35*ss(4);
pos_4=[pl,pb,pw,ph]; pos_5=[0.02*pw,0.8*ph,0.95*pw,0.2*ph]; pos_6=[0.02*pw,0.02*ph,0.95*pw,0.75*ph];

%% Model parameters
a = -2;             % Velocity parameter

nx = 200;           % Number of grid points
dx = 0.05;          % Spatial Resolution

nt = 320;           % Total time steps
dt = 0.025;         % Temporal Resolution

mu = a*(dt/dx);      % Courant number print

%% Allocate arrays
x = zeros(nx,1);           % State
xb = zeros(nx,1);          % X_b Background Initial conditions
xot = zeros(nx,1);         % X_t Truth Initial conditions
xtstore = zeros(nx,nt+1);  % Truth trajectory storage
xastore = zeros(nx,nt+1);  % Analysis trajectory storage
xbstore = zeros(nx,nt+1);  % Background trajectory storage
ystore = zeros(nx,nt+1);   % Observation storage

% Error calc arrays for plots
xerror = zeros(nt+1,1); xrelerror = zeros(nt+1,1); xtnorm = zeros(nt+1,1);
etaerror = zeros(nt,1); etarelerror = zeros(nt,1); etatnorm = zeros(nt,1);

% Physical space index 
xs = zeros(nx,1);
for i=1:nx
xs(i) = (i-1)*dx;
end

%-------------------------------------------------------------------------
%  True error generation
%-------------------------------------------------------------------------
I = eye(nx,nx);
% Generating errors with spatial correlations
Ls = 0.1;                 % Correlation length scale
[C] = SOAR(nx,dx,0.1);    % SOAR Function

% True Model Error Cov matrix
var_qt = 1e-3;            % True model error variance
Qt = var_qt*I; QtInv = inv(Qt); Qtroot = sqrtm(Qt);

qbar = 0;
q = qbar*ones(nx,nt);     % Model error bias
etat = zeros(nx,nt);      % True errors array

% % Constant BIAS Error generation
% etat(1:nx,1:nt) = qbar;

% % Random BIAS Error generation (Random in time, spatially constant model error)
% for i=1:nt
%   etat(:,i) = sqrt(var_qt)*randn(1)*ones(nx,1);
% end

% Error generation (Random spatially correlated model error, with bias)
etat = q + Qtroot*randn(nx,nt);

% Random spatially correlated model error with random bias??

%-------------------------------------------------------------------------
%  2. Assimilation Covariances
%-------------------------------------------------------------------------
var_b = 1e-3; var_o = 1e-3; var_q = 1e-3; % Variances; bckground, obs, model error

B = var_b*I; Binv = inv(B); % Background error covariance matrix
R = var_o*I; Rinv = inv(R); % Observation error covariance matrix
Q = var_q*I; Qinv = inv(Q); % Model error covariance matrix 

%-------------------------------------------------------------------------
%  3. Calculate and store truth trajectory
%-------------------------------------------------------------------------

%% Initial Conditions: Gaussian-shaped curve
xc = (nx*dx)/2; sigma = 1; Ampl = 2; % Mode, Standard dev and Amplitude     
for i=1:nx
    xot(i) = Ampl*exp(-(xs(i)-xc).^2/sigma^2);
end

% %% Truth Trajectory
% xt=xot;
% for n=1:nt
% xtstore(:,n)=xt;                          
% [mx]=adv1dtruth(xt,nx,a,1,dt,dx,etat(:,n)); 
% xt = mx;
%  figure(10), plot(xs,xt,'b'), title('Truth run'), xlabel('Position'), legend('Analysis from t_0 to t_{now}'), ylabel('Value')
%  drawnow
% end
% xtstore(:,nt+1)=xt;

%-------------------------------------------------------------------------
% 3. Set up and store background field with errors
%-------------------------------------------------------------------------

% % Background error perturbation
% bpert = 0; % Background variance noise
% for i=1:nx
%     xb(i) = xot(i) + sqrt(bpert)*randn(1,1);
% end

% Background phase error
bshift = 3;
for i=1:nx
    xb(i) = Ampl*exp(-(xs(i)-xc-bshift).^2/sigma^2);
end

% Background Trajectory
xt=xb;
for n=1:nt
xbstore(:,n)=xt;             
[mx]=adv1d(xt,nx,a,1,dt,dx);
xt=mx;
end
xbstore(:,nt+1)=xt;          % Storage of background trajectory

%---------------------------------------------------------------------
%  4. Take observations of truth model run
%---------------------------------------------------------------------

%% Observation creation parameters
obs_f_t = 1;                    % Obs frequency (time)
obs_f_x = 5;                    % Obs frequency (space)
obsvart = 0;                    % True observation error variance

%% Create observation operator h
h = zeros(nx,nx);

j = 1;
for i = obs_f_x:obs_f_x:nx
    h(i,obs_f_x*j) = 1.;
    j = j + 1;
end

% Observe first and/or final time??

for n=2:nt+1        % This can start from 1 if one wishes to observe first timestep
    
y=zeros(nx,1);      % Setup observation vector

 if rem(n,obs_f_t) == 0   % Obs every obs_f_t timesteps (not 1st time step)    
for i=1:nx
    y(i) = xtstore(i,n) + sqrt(obsvart)*randn(1,1);
end
    ht(n)=1;            
 end

y = h*y; ystore(:,n) = y;

end

obst = sum(ht)                  % Obs total (time)
obsx = sum(sum(h))              % Obs total (space)

%-----------------------------------------------------------------------------%
%(%%) Adjoint & Gradient test                                                 %
% [Forward_Product,Adjoint_Product,Accuracy] = test_adjoint(a,nx,dx,nt,dt)    %
% [phi] = test_grad(xb,h,ht,a,dt,dx,nt,nx,Binv,Rinv,Qinv,ystore);             %
%-----------------------------------------------------------------------------%

% ---------------------------------------------------------------------
% 5. Input minimisation info and perform minimisation
% ---------------------------------------------------------------------

%% Set initial guess
xo = xb; eta = zeros(nx,nt); z = zeros((nt+1)*nx,1);
z(1:nx) = xo; z(nx+1:nx*(nt+1)) = reshape(eta,nx,nt);

%% Minimisation parameters
xstr='Maximum number of iterations'; ystr='Tolerance'; zstr={'30','1d-6'};
min_str=inputdlg({xstr,ystr},'Minimization control',1,zstr); max_iterations=str2num(min_str{1}); tolerance=str2num(min_str{2});

% % Polack-Ribiere Non-Linear Conjugate Gradient
%  [X_a JX dJX it] = minimize_mod_crit(z,'calcJDJ',max_iterations,tolerance, ...
%   xb,h,ht,a,dt,dx,nt,nx,Binv,Rinv,Qinv,ystore);

%% Linear Conjugate gradient (modified matlab pcg giving J on each iteration)
bb = -calcJDJ(zeros(size(z)),xb,h,ht,a,dt,dx,nt,nx,Binv,Rinv,Qinv,ystore);
afun = @(z) calcJDJ(z,xb,h,ht,a,dt,dx,nt,nx,Binv,Rinv,Qinv,ystore) + bb;

[X_a,FLAG,relres,it,dJX,xmat] = minimise_pcg(afun,bb,tolerance,max_iterations,[],[],z);

% %% MinRes minimisation
% bb = -calcJDJ(zeros(size(z)),xb,h,ht,a,dt,dx,nt,nx,Binv,Rinv,Qinv,ystore);
% afun = @(z) calcJDJ(z,xb,h,ht,a,dt,dx,nt,nx,Binv,Rinv,Qinv,ystore) + bb;
% 
% [X_a,FLAG,relres,it,dJX] = minres(afun,bb,tolerance,max_iterations,[],[],z);

% Compute J on each minimisation iterate for convergence plot
for i=1:size(xmat,2)
    zx=xmat(:,i);
    [~,J]=calcJDJ(zx,xb,h,ht,a,dt,dx,nt,nx,Binv,Rinv,Qinv,ystore);
    JX(i)=J;
end

%  X_a contains solution for assimilation window
z = reshape(X_a,nx,nt+1); X_0 = z(:,1); eta = z(:,2:nt+1);

%---------------------------------------------------------------------
% 6. Calculate analysis using forward forecasting model
%    Also calculate error between analysis and truth
%---------------------------------------------------------------------

xerror(1)=norm(xtstore(:,1)-X_0);   % Error norm between truth and analysis
% Now put the solution into forcasting model to get full analysis
xd2=X_0;
for n=1:nt
xastore(:,n)=xd2;
[mx]=adv1d(xd2,a,nx,dx,1,dt);
xd2=mx + eta(:,n);
% figure(10), plot(xs,xd2,'b'), title('Analysis run'), xlabel('Position'), legend('Analysis from t_0 to t_{now}'), ylabel('Value')
% drawnow
%(%&*) While plotting analysis trajectory, calculate 
xerror(n+1)=norm(xtstore(:,n+1)-xd2);
end
xastore(:,nt+1)=xd2;

%---------------------------------------------------------------------
% 7. Plots 
%---------------------------------------------------------------------

%% Convergence plot
h3=figure('Position',pos_6); clf; j = max(size(dJX)); xvals=[0:j-1];
if exist('JX','var') == 1     
j = max(size(JX)); 
subplot(2,1,1)
semilogy(xvals,JX), title('Convergence of cost function'), xlabel('Iteration'), ylabel('Cost function')
subplot(2,1,2)
   grad=zeros(j,1);
for i = 1 : j
  g =  dJX(i,:);
  grad(i)=norm(g);
end
semilogy(xvals,grad), title('Convergence of gradient'), xlabel('Iteration'), ylabel('Norm of gradient')
else
    grad=zeros(j,1);
for i = 1 : j
  g =  dJX(i,:);
  grad(i)=norm(g);
end
semilogy(xvals,grad), title('Convergence of gradient'), xlabel('Iteration'), ylabel('Norm of gradient')
end

%% xt, xa, xb
figure('Position',pos_1)
i = 0; ts = 0;
for n = 1:4
    if n==4 
        i=nt; ts=i*dt;
    else
i=(n-1)*(nt/4)+1; ts = (i-1)*dt;
    end
subplot(2,2,n), hold on, ...
plot(xs,xbstore(:,i),'--g'), title(['t=' num2str(ts)]), ...
plot(xs,xastore(:,i),'r'), plot(xs,xtstore(:,i)), xlabel('x'), ylabel('Value'), ...
end
legend('Background','Analysis','Truth','Location',[0.86 0.85 0.1 0.1])

%% eta^t and eta
figure('Position',pos_3)
i = 0; ts = 0;
for n = 1:4
    if n==4 
i = nt; ts = i*dt;
    else
i = (n-1)*(nt/4) + 1; ts = (i-1)*dt;
    end
subplot(2,2,n), hold on, plot(xs,etat(:,i)), plot(xs,eta(:,i),'r'), xlabel('x'), ...
ylabel('Value'), title(['\eta_{' num2str(i-1) '}']), ...
end
qbar = num2str(mean(mean((eta)))); qbart = num2str(mean(mean((etat)))); qbar = qbar(1:6); qbart = qbart(1:6);
text(-0.5,-0.25,['b=' num2str(qbar) ', b^t=' num2str(qbart)],'sc')
legend('\eta^t \times \Delta t','\eta','Location',[0.86 0.89 0.1 0.1])

%% Error plot (xa vs xt)
time=dt*(0:nt); timeeta=dt*(1:nt);
% eta error calc
for i=1:nt
    etarelerror(i) = 100*(norm(etat(:,i) - eta(:,i)) / norm(etat(:,i)));
end

% state error calc
for i=1:nt+1
    xrelerror(i) = 100*(xerror(i) / norm(xtstore(:,i)));
end

figure('Position',pos_2), clf, plot(time,xrelerror), hold on, ...
    plot(timeeta,etarelerror,'r'), title('Error plot; Truth vs Analysis'), ...
    xlabel('Time (s)'), ylabel('Error Norm (%)'), ...
    legend('State', '\eta', 'Location', [0.89 0.88 0.1 0.1])
drawnow

% %% Covariance Matrix Q
% figure('Position',pos_5), contour(Q), colorbar, set(gca,'YDir','reverse'), xlabel('Model gridpoint'), ylabel('Model gridpoint'), title('Covariance Matrix Q used in Assimilation')

%% Options print out
min_str{3} = num2str(nx*dx);
s1=num2str(dt); min_str{5} = s1; s2=num2str(dx); min_str{6} = s2;
min_str{7} = num2str(obs_f_t); min_str{8} = num2str(obs_f_x); min_str{9} = num2str(obsvart); 
min_str{10} = num2str(bshift); % min_str{10} = num2str(bpert);
min_str{12} = num2str(numel(JX)); min_str{13} = num2str(a);
min_str{14} = num2str(var_o); min_str{15} = num2str(var_b); min_str{16} = num2str(var_q); 
min_str{18} = num2str(obsx); min_str{19} = num2str(nx); min_str{20} = num2str(obst); min_str{21} = num2str(nt); 
s3=num2str(mu);

if numel(s3)<5
min_str{17}= s3;
else
     min_str{17}= s3(1:5);
end

h4=figure('Position',pos_4); clf;
text1 = {['List of options chosen']};
text2 = {['Length of assimilation window: ', min_str{3} 's']};
text3 = {['dt: ' min_str{5}, '   dx: ' min_str{6}, '   a: ' min_str{13}, '   C: ' min_str{17}]};
text4 = {['Observations; ', 'Freq (t): ' min_str{7}, '  Total (t): ' min_str{20} ' / ' min_str{21}, '  Freq (spatial): ' min_str{8}, '  Total (spatial): ' min_str{18} ' / ' min_str{19}]};
text5 = {['Observations true noise variance = ', min_str{9}]};
text6 = {['Background true noise variance = ', min_str{10}]};
%text6={['Background phase shift (in truth) = ', min_str{10}]};
text7 = {['Background Cov:  ' min_str{15}, '   Observations Cov:  ' min_str{14}, '   Model Error Cov:  ' min_str{16}]};
text8 = {['Maximum iterations: ', min_str{1}]};
text9 = {['Number of actual iterations: ', min_str{12}]}; 
text10 = {['Tolerance: ', min_str{2}]};

str1 = [text2;text3;text4;text5;text6;text7;text8;text9;text10];
uicontrol('Style','text','Position',pos_5,'String',text1,...
 'FontSize',14,'FontWeight','bold')
uicontrol('Style','text','Position',pos_6,'String',str1,...
 'FontSize',12,'HorizontalAlignment','left')