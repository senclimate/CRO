function [T_std,h_std]=RO_std_anal(par)
%% Solve Analytical Standard Deviation of the Linear-White_Additive RO
% White Noise RO (ntype==1)
% dT/dt=R*T+F1*h+(sigma_T)*(w_T)
% dh/dt=-F2*T-epsilon*h+(sigma_h)*(w_h)
% where w_T and w_h are white noise with zero mean and unit variance

%% Setting Parameter
% post-processing par
for i=1:length(par)
    par_sel(i)=par{i}(1); % only select annual mean value
end 

% linear parameter
R=par_sel(1);
F1=par_sel(2);
F2=par_sel(3);
epsilon=par_sel(4);

% noise parameter
sigma_T=par_sel(9);
sigma_h=par_sel(10);

%% Calculate Standard Deviation
T_std=sqrt(((F1*F2-epsilon*R+epsilon^2)*(sigma_T^2)+(F1^2)*(sigma_h^2))/(2*(-R+epsilon)*(F1*F2-R*epsilon)));
h_std=sqrt(((F2^2)*sigma_T^2+(F1*F2-epsilon*R+R^2)*sigma_h^2)/(2*(-R+epsilon)*(F1*F2-R*epsilon)));

end 