function [T_anal_out,h_anal_out,noise_out]=RO_solver_analytic(par,IC,N,NE,dt,saveat,savemethod,noise_custom)
%% Solve Analytical Solution for Linear White RO 
% White Noise RO (ntype==1)
% dT/dt=R*T+F1*h+(sigma_T)*(xi_T)
% dh/dt=-F2*T-epsilon*h+(sigma_h)*(xi_h)
% where xi_T and xi_h are white noise with zero mean and unit variance

%% Setting Default Option
if nargin<5||isempty(dt)
    dt=0.1; 
end
if nargin<6||isempty(saveat)
    saveat=1.0;  
end
if nargin<7||isempty(savemethod)
    savemethod="sampling";
end
if nargin<8||isempty(noise_custom)
    noise_custom=[];  
end

%% Setting Time
NT=round(N/dt);

%% Setting Noise
if isempty(noise_custom)
    noise_custom=randn(NT-1,4,NE);
else
    if isscalar(noise_custom)
        rng(noise_custom)
        noise_custom=randn(NT-1,4);
        noise_custom=repmat(noise_custom,[1 1 NE]); % note: must keep consistency with noise generator in RO_solver        
    end 
end 

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

%% Setting Values 
To=IC(1);
ho=IC(2);

t=0:dt:(NT-1)*dt;
t=transpose(t); 

gr=(R-epsilon)/2;
w=sqrt(4*F1*F2-(R+epsilon)^2)/2;

%% Calculate Deterministic Part
T_anal_det=exp(gr*t).*(To*cos(w*t)+To*(R+epsilon)/(2*w)*sin(w*t)+ho*(F1/w)*sin(w*t));
h_anal_det=exp(gr*t).*(ho*cos(w*t)-ho*(R+epsilon)/(2*w)*sin(w*t)-To*(F2/w)*sin(w*t));

%% Calculate Stochastic Part
for j=1:NE
w_T=sigma_T*squeeze(noise_custom(:,1,j))/sqrt(dt);
w_h=sigma_h*squeeze(noise_custom(:,2,j))/sqrt(dt);

T_anal_sto=zeros(length(t),1);
h_anal_sto=zeros(length(t),1);

for i=1:length(t)
    if i==1
        T_anal_sto(1,1)=0.0;
        h_anal_sto(1,1)=0.0;
    else
        tau=t(1:i-1);
        T_anal_sto(i,1)=sum(dt*exp(gr*(t(i)-tau)).*(w_T(1:i-1).*cos(w*(t(i)-tau))+w_T(1:i-1).*(R+epsilon)/(2*w).*sin(w*(t(i)-tau))+w_h(1:i-1).*(F1/w).*sin(w*(t(i)-tau))),1);
        h_anal_sto(i,1)=sum(dt*exp(gr*(t(i)-tau)).*(w_h(1:i-1).*cos(w*(t(i)-tau))-w_h(1:i-1).*(R+epsilon)/(2*w).*sin(w*(t(i)-tau))-w_T(1:i-1).*(F2/w).*sin(w*(t(i)-tau))),1);
    end
end 

%% Total Solution 
T_anal=T_anal_det+T_anal_sto;
h_anal=h_anal_det+h_anal_sto;

T_anal_out(:,j)=T_anal;
h_anal_out(:,j)=h_anal;
end 

%% Select Time Step % Save Noise Output
% note: noise_out does not follow this savemethod option 
if savemethod=="sampling"
    T_anal_out=T_anal_out(1:round(saveat/dt):end,:);
    h_anal_out=h_anal_out(1:round(saveat/dt):end,:);
elseif savemethod=="mean"
    T_out_sel=T_anal_out(0.5*round(1/dt):end-0.5*round(1/dt)-1,:);
    h_out_sel=T_anal_out(0.5*round(1/dt):end-0.5*round(1/dt)-1,:);  
    T_out_sel=reshape(T_out_sel,round(saveat/dt),[],size(T_out_sel,2));
    T_out_sel=squeeze(mean(T_out_sel,1));
    h_out_sel=reshape(h_out_sel,round(saveat/dt),[],size(h_out_sel,2));
    h_out_sel=squeeze(mean(h_out_sel,1));
    if NE>1
    T_anal_out=[T_anal_out(1,:);T_out_sel];
    h_anal_out=[h_anal_out(1,:);h_out_sel];
    elseif NE==1
    T_anal_out=[T_anal_out(1),T_out_sel];
    h_anal_out=[h_anal_out(1),h_out_sel];
    T_anal_out=transpose(T_anal_out);
    h_anal_out=transpose(h_anal_out);    
    end 
end
noise_out=noise_custom; 

end 