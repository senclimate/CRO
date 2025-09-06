function [f_T,g_T,f_h,g_h,f_xi_T,g_xi_T,f_xi_h,g_xi_h]=RO_tendency(par,T,h,xi_T,xi_h,EF)
%% setting parameter
% linear parameter
R=par(1);
F1=par(2);
F2=par(3);
epsilon=par(4);

% nonlinear parameter
b_T=par(5);
c_T=par(6);
d_T=par(7);
b_h=par(8);

% noise parameter
sigma_T=par(9);
sigma_h=par(10);
B=par(11);
m_T=par(12);
m_h=par(13);

% noise option parameter
n_T=par(14);
n_h=par(15); 
n_g=par(16); 

%% setting external forcing
E_T=EF(1);
E_h=EF(2);

%% setup noise type
if n_g==0
    g_T_factor=B*T;
elseif n_g==1
    h_T_factor=T; 
    h_T_factor(T<=0)=0; % faster than heaviside(T)
    g_T_factor=B*h_T_factor; 
elseif n_g==2
    g_T_factor=0;
end 

%% calculate tendency
% dT/dt
if n_T==0 % red noise
    f_T=R*T+F1*h+b_T*(T^2)-c_T*(T^3)+d_T*T*h+sigma_T*(1+g_T_factor)*xi_T+E_T;
    g_T=0.0;    
elseif n_T==1 % white noise
    f_T=R*T+F1*h+b_T*(T^2)-c_T*(T^3)+d_T*T*h+E_T;
    g_T=sigma_T*(1+g_T_factor);
else
    disp("Error: Please specify valid n_T value (0 or 1)")
    f_T=NaN;
    g_T=NaN;
end

% dh/dt
if n_h==0 % red noise
    f_h=-F2*T-epsilon*h-b_h*(T^2)+sigma_h*xi_h+E_h;
    g_h=0.0;    
elseif n_h==1 % white noise
    f_h=-F2*T-epsilon*h-b_h*(T^2)+E_h;
    g_h=sigma_h;
else
    disp("Error: Please specify valid n_h value (0 or 1)")
    f_h=NaN;
    g_h=NaN;    
end

% d(xi_T)/dt
if n_T==0 % red noise
    f_xi_T=-m_T*xi_T;
    g_xi_T=sqrt(2*m_T); 
elseif n_T==1 % white noise (this part doens't use for calculation)
    f_xi_T=NaN;
    g_xi_T=NaN;
end

% d(xi_h)/dt
if n_h==0 % red noise
    f_xi_h=-m_h*xi_h;
    g_xi_h=sqrt(2*m_h); 
elseif n_h==1 % white noise (this part doens't use for calculation)
    f_xi_h=NaN;
    g_xi_h=NaN;    
end

end