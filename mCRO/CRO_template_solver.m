clc
clear

%%-------------------------------------------------------------------------
% This template allows users to easily configure options for solving the RO model
%%-------------------------------------------------------------------------
% RO model Equations
%--------------------------------------------------------------------------
% dT/dt=R*T+F1*h+b_T*(T^2)-c_T*(T^3)+d_T*(T*h)+sigma_T*(1+B*G(T)*T)*noise_T
% dh/dt=-F2*T-epsilon*h-b_h*(T^2)+sigma_h*noise_h
%--------------------------------------------------------------------------
% If noise is white: 
% noise_T=w_T (n_T=1)
% noise_h=w_h (n_h=1)
%--------------------------------------------------------------------------
% If noise is red: 
% d(noise_T)/dt=-m_T*T+sqrt(m_T)*w_T (n_T=0)
% d(noise_h)/dt=-m_h*h+sqrt(m_h)*w_h (n_h=0)
%--------------------------------------------------------------------------
% If T noise type is: 
% Additive: G(T)=1 (n_g=2)
% Multiplicative-Linear: G(T)=B*T (n_g=0) 
% Multiplicative-Heaviside: G(T)=B*H(T)*T (n_g=1)

%% Setting RO Parameter
% Choose one of the following to set 'par' (Method 1, 2 and 3)  

%% Method 1: Use Parmaeter Libarary Set (*replace with the option the user wants to run)
data_name="ORAS5"; % data type of fitted parameters (see CRO manual for full list) 
RO_type="Linear-White-Additive"; % RO type of fitted parameters (see CRO manual for full list) 
par=par_load(data_name,RO_type); 

%% Method 2: Use User-Specified Parameter Set (*replace with the parameter the user wants to run)
% Linear Parameter 
R=-0.05; % -0.05 dT/dt=R*T
F1=0.02; % dT/dt=F1*h
F2=0.9; % dh/dt=-F2*T
epsilon=0.03; % dh/dt=-epsilon*h

% Nonlinear Parameter
b_T=0.0; % dT/dt=(b_T)*(T^2)
c_T=0.0; % dT/dt=-(c_T)*(T^3)
d_T=0.0; % dT/dt=(d_T)*(T*h)
b_h=0.0; % dh/dt=-(b_h)*(T^2)

% Noise Parameter 
sigma_T=0.2; % dT/dt=(sigma_T)*(N_T) 
sigma_h=1.2; % dh/dt=(sigma_h)*(N_h) 
B=0.0; % dT/dt=(sigma_T)*(1+B*T)*(N_T) or dT/dt=(sigma_T)*(1+B*H(T)*T)*(N_T)
m_T=1.0; % d(xi_T)/dt=-m*T*(xi_T);
m_h=1.0; % d(xi_h)/dt=-m*h*(xi_h);

% Noise Option Parameter 
n_T=1; % noise type for T (0: red noise, 1: white noise)
n_h=1; % noise type for h (0: red noise, 1: white noise)
n_g=0; % multiplicative noise type for T (0: linear, 1: Heaviside linear, 2: omit this option) (note: only valid when B is not zero) 

% Save into cell array 
par={R;F1;F2;epsilon;b_T;c_T;d_T;b_h;sigma_T;sigma_h;B;m_T;m_h;n_T;n_h;n_g};

% Example: Impose Seasonal Varition in R
par{1}=[R 2*R pi 0.5*R 0]; % X, Xa, φa, Xsa, φsa (X=X+Xa*sin(wt+φa)+Xsa*sin(2wt+φsa), where w=2π/12). Alternatively, par{1}=[R 2*R pi]; 

%% Method 3: Use returned by 'RO_fitting' function 
% See 'CRO_template_fitting.m' to see how to obtain the 'par' using 'RO_fitting' using your own T and h data for a user-selected RO type. 

%% Setting RO Solver Options (*replace with the option the user wants to run)
IC=[1.0;0.0]; % initial condition of T (K) and h (m) 
N=12*(10); % simulation time (month) 
NE=100; % number of ensemble 
NM="EH"; % numerical method of RO solver ("EH" and "EM") (optional input; default="EH") 
dt=0.1; % numerical time step (month) (optional input; default=0.1) 
saveat=1.0; % save interval of simulated T and h (month) (optional input; default=1.0)
savemethod="sampling"; % save method ("sampling" and "mean") (optional input, default="sampling")
EF={0.0;0.0}; % external forcing (optional input; deafult={0.0 0.0})
noise_custom=[]; % customized noise for solving RO (optional input; default=[]) 

%% Solve RO and Plot
% Solve
[T,h]=RO_solver(par,IC,N,NE,NM,dt,saveat,savemethod,EF,noise_custom); % output: T ([time x ensemble]), h ([time x ensemble])

% Plot
time=0:saveat:(round(N/saveat)-1)*saveat;

% Plot First Ensemble T
figure;
plot(time,T(:,1))

% Plot Ensemble-mean T 
figure;
plot(time,mean(T,2))
