clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Weak-Constraint 4D-Var with model error forcing applied to 
%  the 1D Advection model
%
%  Written by Adam El-Said
%
%  List of main variables
%
%% Assimilation variables
%    xb:         Background term
%    h:          Spatial component of observation operator
%    ht:         Temporal component of observation operator
%    B:          Background weighting matrix (inversed)
%    R:          Observation weighting matrix (inversed)
%    Q:          Model Error weighting matrix (inversed)
%    [z]:        Current iterate
%    [ystore]:   Observation values
%
%% Model variables
%    L:          Physical domain size
%    a:          Velocity in advection equation
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

% Close any open figure windows
g=gcf;
for i=1:g
  close(i)
end

%% Graph position settings
ss=get(0,'ScreenSize');
fig_width=0.55*ss(3);
fig_height=0.4*ss(4);
pos_1=[0.02*ss(3),0.5*ss(4),fig_width,fig_height];
pos_2=[0.02*ss(3),0.05*ss(4),fig_width,fig_height];
pos_3=[0.6*ss(3),0.45*ss(4),0.4*ss(3),0.45*ss(4)];
pl=0.6*ss(3);
pb=0.05*ss(4);
pw=0.3*ss(3);
ph=0.35*ss(4);
pos_4=[pl,pb,pw,ph];
pos_5=[0.02*pw,0.8*ph,0.95*pw,0.2*ph];
pos_6=[0.02*pw,0.02*ph,0.95*pw,0.75*ph];

%% Model parameters
a = -2;                  % Velocity parameter

L = 10;                  % Domain size
dx = 0.05;               % Spatial Resolution
nx = L/dx;               % Number of grid points

tassim = 8;              % Assimilation window time length
dt = 0.025;              % Temporal Resolution
nt = tassim/dt;          % Total time steps

C = a*(dt/dx)            % Courant number print

%% Allocate arrays
x = zeros(nx,1);         % State
xb = zeros(nx,1);        % X_b Background Initial conditions
xot = zeros(nx,1);       % X_t Truth Initial conditions
xtstore = zeros(nx,nt);  % Truth trajectory storage
xastore = zeros(nx,nt);  % Analysis trajectory storage
xbstore = zeros(nx,nt);  % Background trajectory storage
xerror = zeros(nt,1);    % Error calculation between analysis and truth

% Physical space index
xs = zeros(nx,1);
for i=1:nx
xs(i) = i*dx;
end

%-------------------------------------------------------------------------
%  True error generation
%-------------------------------------------------------------------------

% Generating errors with spatial correlations
Ls = 0.1;                         % Correlation length scale
[C] = SOAR(nx,Ls,L);              % SOAR Function
% C = eye(nx);                    % Spatially uncorrelated
Cinv = inv(C);
sd_qt = 1;                        % True model error sd
Sigma_qt = sd_qt*eye(nx);         % Diagonal matrix, model error sd's

Qt = Sigma_qt*C*Sigma_qt;         % Auto-Cov matrix now ready
QtInv = inv(Qt);                  % Inverse
Qtroot = sqrtm(Qt);               % Square root

b = 0;
mut = b*ones(nx,1);               % True model error mean/bias
etat = zeros(nx,nt-1);            % True errors array

% % Error generation (Constant model error)
% etat = zeros(nx,nt-1);
% etat(:,:) = 5e-2;

% % Error generation (Random in time, spatially constant model error)
% for i=1:nt-1
%     etat(:,i) = mut + sd_qt*randn(1)*ones(nx,1);
% end

% Error generation (Random spatially correlated model error, with bias)
% **(For Spatially uncorrelated model error, set C = I)**
for i=1:nt-1
etat(:,i)= mut + Qtroot*randn(nx,1);
end

etay = etat*dt;
%-------------------------------------------------------------------------
%  2. Assimilation Covariances
%-------------------------------------------------------------------------

var_b = 1e-3;                 % Background variance and standard dev
var_o = 1e-3;                 % Observations variance
var_q = 1e-3;                 % Model Error variance

sd_b = sqrt(var_b);           % Background standard deviation                 
sd_o = sqrt(var_o);           % Observations standard deviation
sd_q = sqrt(var_q);           % Model Error standard deviation

Sigma_b = eye(nx)*(1/sd_b);   % Diagonal Background sd error matrix
Sigma_o = eye(nx)*(1/sd_o);   % Diagonal Observation sd error matrix
Sigma_q = eye(nx)*(1/sd_q);   % Diagonal Model sd error matrix

B = Sigma_b*eye(nx)*Sigma_b;  % Background error covariance matrix (inversed)
R = Sigma_o*eye(nx)*Sigma_o;  % Observation error covariance matrix (inversed)
Q = Sigma_q*Cinv*Sigma_q;     % Model error covariance matrix (inversed)

%-------------------------------------------------------------------------
%  3. Calculate and store truth trajectory
%-------------------------------------------------------------------------

%% Initial Conditions are a Gaussian-shaped curve
xc = L/2;                     % Domain midpoint (mean/mode)
sigma = 1;                    % Variance
Ampl = 2;                     % Amplitude
for i=1:nx
    xot(i) = Ampl*exp(-(xs(i)-xc).^2/sigma^2);
end

%% Truth Trajectory
xt=xot;
for n=1:nt-1
xtstore(:,n)=xt;                           % Storage of trajectory
[mx]=adv1dtruth(xt,L,a,1,dt,dx,etat(:,n)); % 1 time step
xt=mx;
end
xtstore(:,nt)=xt;

%-------------------------------------------------------------------------
% 3. Set up and store background field with errors
%-------------------------------------------------------------------------

% Background phase error
bshift = 0;
for i=1:nx
    xb(i) = Ampl*exp(-(xs(i)-xc-bshift).^2/sigma^2);
end

% % Background error perturbation
% bpert = 0.01; % Background variance noise
% for i=1:nx
%     xb(i) = xot(i) + sqrt(bpert)*randn(1,1);
% end

% Background Trajectory
xt=xb;
for n=1:nt-1
xbstore(:,n)=xt;                      % Storage of trajectory
[mx]=adv1d(xt,0*eye(nx),L,a,1,dt,dx); % 1 time step
xt=mx;
end
xbstore(:,nt)=xt;

%---------------------------------------------------------------------
%  4. Take observations of truth model run
%---------------------------------------------------------------------

%% Observation creation parameters
ystore = zeros(nx,nt);      % Observation storage

obs_f_t = 5;                % Obs frequency  (time)
obs_f_x = 7;                % Obs frequency  (space)

obst = floor(nt/obs_f_t);   % Obs total (time)
obsx = floor(nx/obs_f_x);   % Obs total (space)
obspert = 0;                % Error variance of obs from the truth

%% Create observation operator h
h = zeros(nx,nx);

j=1;
for i=obs_f_x:obs_f_x:nx
    h(i,obs_f_x*j)=1.;
    j=j+1;
end

ht = zeros(nt,1); % Observation time vector. 
                  % Put simply =1 when an obs is taken, 0 otherwise. 
xd=xot;
for n=1:nt

y=zeros(nx,1);            % Setup observation vector

 if rem(n,obs_f_t) == 0   % Obs every obs_f_t timesteps (not 1st time step)
     for i=1:nx
    y(i) = xtstore(i,n) + sqrt(obspert)*randn(1,1);
     end
     ht(n)=1;             %(*) =1 when an observation was taken
 end

y = h*y;          % This operation means y is observed every obs_f_x grid point
ystore(:,n) = y;  % Store all y observation vectors

end

% ---------------------------------------------------------------------
% (%%) Gradient test - To run this, put a break point after the above is
%                      complete and run function test_adv1d_grad manually.
% ---------------------------------------------------------------------

% [phi] = test_adv1d_grad(xb,h,ht,L,a,dt,dx,nt,nx,B,R,Q,ystore);

% ---------------------------------------------------------------------
% 5. Input minimization info and perform minimization
% ---------------------------------------------------------------------

%% Set initial values of J variables
xo = xb;                 % x^a_0 = x^b_0
eta = zeros(nx,nt-1);    % Error increment eta = 0 initially
z = zeros(nt*nx,1);      % x0 and all eta_i's stored here for gradient solving procedure              

% z=[xo,eta(1),...,eta(nx)]
% xo and all eta's are arranged ready for minimisation

k=1;
      for i=1:nx
          z(k) = xo(i);
          k=k+1;
      end

k=nx+1;
    for n=1:nt-1
      for i=1:nx
          z(k) = eta(i,n);
          k=k+1;
      end
    end

xstr='Maximum number of iterations';
ystr='Tolerance';
zstr={'30','1d-5'};
min_str=inputdlg({xstr,ystr},'Minimization control',1,zstr);
max_iterations=str2num(min_str{1});
tolerance=str2num(min_str{2});

[X_a JX dJX it] = minimize_mod_crit(z,'calcJDJ',max_iterations,tolerance, ...
xb,h,ht,L,a,dt,dx,nt,nx,B,R,Q,ystore);

%  X_a will have the contain the solution for the full assimilation window(
% (X_0 and all ETA's). So we need to redistribute them in the correct form.

X_0 = X_a(1:nx);
 n=1; i=1;
for k=(nx+1):(nx*nt)
 eta(i,n) = X_a(k);
 i = i+1;
  if rem(k,nx) == 0
    n=n+1;
    i=1;
  end
end

%---------------------------------------------------------------------
% 6. Calculate analysis using forward forecasting model
%    Also calculate error between analysis and truth
%---------------------------------------------------------------------

xerror(1)=norm(xtstore(:,1)-X_0);   % Error norm between truth and analysis
% Now put the solution into forcasting model to get full analysis
xd2=X_0;
for n=1:nt-1
xastore(:,n)=xd2;
[X_A]=adv1d(xd2,eta(:,n),L,a,1,dt,dx);
xd2=X_A;
%(%&*) While plotting analysis trajectory, calculate 
xerror(n+1)=norm(xtstore(:,n+1)-xd2);
end
xastore(:,nt)=xd2;

% eta error calc
for i=1:nt-1
    etaerror(i) = norm( etay(:,i) - eta(:,i));
    etanorm(i) = norm(etat(:,i));
    etarelerror(i) = 100*(etaerror(i) / etanorm(i));
end

% state error calc
for i=1:nt
    xtnorm(i) = norm(xtstore(:,i));
    xrelerror(i) = 100*(xerror(i) / xtnorm(i));
end

%---------------------------------------------------------------------
% 7. Plots 
%---------------------------------------------------------------------

%% xt, xa, xb
figure('Position',pos_1)
i = 0; ts = 0;
for n = 1:4
    if n==4 
        i=nt-1;
        ts=(i+1)*dt;
    else
i=(n-1)*(nt/4)+1;
ts = (i-1)*dt;
    end
 subplot(2,2,n), hold on, ...
    plot(xs,xbstore(:,i),'--g'), title(['t=' num2str(ts)]), plot(xs,xastore(:,i),'r'), plot(xs,xtstore(:,i)), ...
     xlabel('Position (x)'), ylabel('Flux'), ...
end
legend('Background','Analysis','Truth','Location',[0.86 0.85 0.1 0.1])

%% eta truth and eta
figure('Position',pos_3)
i = 0; ts = 0;
for n = 1:4
    
    if n==4 
        i=nt-1;
        ts=(i+1)*dt;
    else
i=(n-1)*(nt/4)+1;
ts = (i-1)*dt;
    end

 subplot(2,2,n), hold on, ...
     title(['t=' num2str(ts)]), plot(xs,etay(:,i)), plot(xs,eta(:,i),'r'), ...
     xlabel('Position (x)'), ylabel('Flux'), ...
end
legend('\eta^t \times \Delta t','\eta','Location',[0.86 0.89 0.1 0.1])

%% Convergence plot
h3=figure('Position',pos_6); clf; j = max(size(JX)); xvals=[0:j-1];
subplot(2,1,1)
semilogy(xvals,JX), title('Convergence of cost function'), xlabel('Iteration'), ylabel('Cost function')
subplot(2,1,2)
grad=zeros(j,1);
for i = 1 : j
  g =  dJX(:,i);
  grad(i)=norm(g);
end
semilogy(xvals,grad), title('Convergence of gradient'), xlabel('Iteration'), ylabel('Norm of gradient')

%% Error plot (xa,etaa vs xt,etat)
time=dt*(1:nt);
timeeta=dt*(1:nt-1);
figure('Position',pos_2), clf, plot(time,xrelerror), hold on, ...
    plot(timeeta,etarelerror,'r'), title('Error plot; Truth vs Analysis'), ...
    xlabel('Time (s)'), ylabel('Error Norm (%)'), ...
    legend('State', '\eta', 'Location', [0.89 0.88 0.1 0.1])
drawnow

%% Covariance Matrix Q
QQ = inv(Q); % This is the inverse of the inverse (so its just Q)
figure('Position',pos_5), contourf(QQ), colorbar, set(gca,'YDir','reverse'), ...
    xlabel('Model gridpoint'), ylabel('Model gridpoint'), ...
    title('Covariance Matrix Q used in Assimilation')

%% Options print out
min_str{3} = num2str(tassim);  % min_str{4} = num2str(fcstep);
s1=num2str(dt); min_str{5} = s1; s2=num2str(dx); min_str{6} = s2;
min_str{7} = num2str(obs_f_t); min_str{8} = num2str(obs_f_x); min_str{9} = num2str(obspert); 
min_str{10} = num2str(bshift); % min_str{10} = num2str(bpert); 
min_str{12} = num2str(numel(JX)); min_str{13} = num2str(a);
min_str{14} = num2str(var_o); min_str{15} = num2str(var_b); min_str{16} = num2str(var_q); 
min_str{18} = num2str(obsx); min_str{19} = num2str(nx); min_str{20} = num2str(obst); min_str{21} = num2str(nt); 
s3=num2str(C);

if numel(s3)<5
min_str{17}= s3;
else
     min_str{17}= s3(1:5);
end

h4=figure('Position',pos_4);
clf;
text1={['List of options chosen']};
text2={['Length of assimilation window: ', min_str{3} 's']};
text3={['dt: ' min_str{5}, '   dx: ' min_str{6}, '   a: ' min_str{13}, '   C: ' min_str{17}]};
text4={['Observations; ', 'Freq (t): ' min_str{7}, '  Total (t): ' min_str{20} ' / ' min_str{21}, '  Freq (spatial): ' min_str{8}, '  Total (spatial): ' min_str{18} ' / ' min_str{19}]};
text5={['Observations noise variance (in truth) = ', min_str{9}]};
%text6={['Background noise variance (in truth) = ', min_str{10}]};
text6={['Background phase shift (in truth) = ', min_str{10}]};
text7={['Background Cov:  ' min_str{15}, '   Observations Cov:  ' min_str{14}, '   Model Error Cov:  ' min_str{16}]};
text8={['Maximum iterations: ', min_str{1}]};
text9={['Number of actual iterations: ', min_str{12}]}; 
text10={['Tolerance: ', min_str{2}]};

str1=[text2;text3;text4;text5;text6;text7;text8;text9;text10];
uicontrol('Style','text','Position',pos_5,'String',text1,...
 'FontSize',14,'FontWeight','bold')
uicontrol('Style','text','Position',pos_6,'String',str1,...
 'FontSize',12,'HorizontalAlignment','left')