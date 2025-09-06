clc
clear

%%-------------------------------------------------------------------------
% This template allows users to easily configure options for fitting the RO model
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

%% Load Data to Fit 
% Load T and h (replace with the data the user wants to fit)
fpath=pwd+"/Data/"+"XRO_indices_oras5.nc"; % file load path 
T=ncread(fpath,'Nino34'); % load T (Nino 3.4 index)
h=ncread(fpath,'WWV'); % load h (equatorial mean thermocline depth) 
ncdisp(fpath)

% Plot T 
time=1958:(1/12):1958+46-1/12; % 1958 (Jan) - 2020 (Dec) 

figure;
plot(time,T)

%% Fitting Observation T and h 
% Set RO type (change to the type the user wants to fit)
par_option_T=struct("R",1,"F1",1,"b_T",0,"c_T",0,"d_T",0);
par_option_h=struct("F2",1,"epsilon",1,"b_h",0);
par_option_noise=struct("T","white","h","white","T_type","additive");

% Parameter fitting 
par=RO_fitting(T,h,par_option_T,par_option_h,par_option_noise,"LR-F");


