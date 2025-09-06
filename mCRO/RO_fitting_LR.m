function par=RO_fitting_LR(T,h,par_option_T,par_option_h,par_option_noise,dt,tend_option,fitting_option_B,fitting_option_red)
warning('off')

% Convert to Vector  
T=reshape(T,[length(T),1]);
h=reshape(h,[length(h),1]);

% Convert to Array 
par_option_T=reshape(par_option_T,[1,length(par_option_T)]);
par_option_h=reshape(par_option_h,[1,length(par_option_h)]);
par_option_noise=reshape(par_option_noise,[1,length(par_option_noise)]);

n_T=par_option_noise(1,1); % 0: red noise, 1: white noise
n_h=par_option_noise(1,2); % 0: red noise, 1: white noise
n_g=par_option_noise(1,3); % 0: multi-linear, 1: multi-Heaviside linear, 2: skip this option

% Setting Option
par_option_T_onoff=par_option_T;
par_option_h_onoff=par_option_h;
par_option_T_onoff(find(par_option_T_onoff>=1))=1;
par_option_h_onoff(find(par_option_h_onoff>=1))=1;

par_option_T_season=find(par_option_T>1); % Index for Seasonal Parameter 
par_option_h_season=find(par_option_h>1); % Index for Seasonal Parameter 

% Calculate Differential and Select T and h for fitting 
if tend_option=="F"
Tdiff=diff(T)/dt;
hdiff=diff(h)/dt;
Ts=T(1:end-1);
hs=h(1:end-1);
elseif tend_option=="C"
Tdiff=(T(3:end)-T(1:end-2))/(2*dt);
hdiff=(h(3:end)-h(1:end-2))/(2*dt);
Ts=T(2:end-1);
hs=h(2:end-1);    
end 

% Setting Array 
Tv=[Ts,hs,Ts.^2,-Ts.^3,Ts.*hs].*par_option_T_onoff;
hv=[-Ts,-hs,-Ts.^2].*par_option_h_onoff;

% Setting Seasonal Array 
w=(2*pi)/12;
t=0:dt:(length(T)-1)*dt;
t=transpose(t);

if tend_option=="F"
    t=t(1:end-1);
elseif tend_option=="C"
    t=t(2:end-1);
end 

if tend_option=="F"
    t_shift=t+0.5*dt; % for seasonal parameter estimation
elseif tend_option=="C"
    t_shift=t; % for seasonal parameter estimation
end

if ~isempty(par_option_T_season)
    for i=1:length(par_option_T_season)
        if par_option_T(par_option_T_season(i))==3 % Annual
            Tv_add=[Tv(:,par_option_T_season(i)).*sin(w*t_shift),Tv(:,par_option_T_season(i)).*cos(w*t_shift)];
            Tv=[Tv,Tv_add];
        elseif par_option_T(par_option_T_season(i))==5 % Annual + Semi-Annual
            Tv_add=[Tv(:,par_option_T_season(i)).*sin(w*t_shift),Tv(:,par_option_T_season(i)).*cos(w*t_shift),Tv(:,par_option_T_season(i)).*sin(2*w*t_shift),Tv(:,par_option_T_season(i)).*cos(2*w*t_shift)];
            Tv=[Tv,Tv_add];
        end
    end
end

if ~isempty(par_option_h_season)
    for i=1:length(par_option_h_season)
        if par_option_h(par_option_h_season(i))==3 % Annual
            hv_add=[hv(:,par_option_h_season(i)).*sin(w*t_shift),hv(:,par_option_h_season(i)).*cos(w*t_shift)];
            hv=[hv,hv_add];
        elseif par_option_h(par_option_h_season(i))==5 % Annual + Semi-Annual
            hv_add=[hv(:,par_option_h_season(i)).*sin(w*t_shift),hv(:,par_option_h_season(i)).*cos(w*t_shift),hv(:,par_option_h_season(i)).*sin(2*w*t_shift),hv(:,par_option_h_season(i)).*cos(2*w*t_shift)];
            hv=[hv,hv_add];
        end
    end
end

%% Determinsitic Parameters
% Perform Fitting for Deterministic Part 
[par_T,~,res_T]=regress_std(double(Tdiff),double(Tv));
[par_h,~,res_h]=regress_std(double(hdiff),double(hv));

par_T=num2cell(par_T);
par_h=num2cell(par_h);

% Perform Post-Processing for the Seasonal Part
if ~isempty(par_option_T_season)
    ind_base=length(par_option_T);
    for i=1:length(par_option_T_season)
        if par_option_T(par_option_T_season(i))==3 % Annual
            ind_sin=ind_base+1;
            ind_cos=ind_base+2;
            A=sqrt(par_T{ind_sin}^2+par_T{ind_cos}^2);
            phi=mod(atan2(par_T{ind_cos},par_T{ind_sin}),2*pi);
            if sign(par_T{par_option_T_season(i)}(1))~=sign(A)
                A=-A;
                phi=wrapToPi(phi-pi);
            end 
            par_T{par_option_T_season(i)}(2:3)=[A,phi];
            ind_base=ind_cos;
        elseif par_option_T(par_option_T_season(i))==5 % Annual + Semi-Annual
            ind_sin_a=ind_base+1;
            ind_cos_a=ind_base+2;
            ind_sin_sa=ind_base+3;
            ind_cos_sa=ind_base+4;        
            A_a=sqrt(par_T{ind_sin_a}^2+par_T{ind_cos_a}^2);
            phi_a=mod(atan2(par_T{ind_cos_a},par_T{ind_sin_a}),2*pi);
            A_sa=sqrt(par_T{ind_sin_sa}^2+par_T{ind_cos_sa}^2);
            phi_sa=mod(atan2(par_T{ind_cos_sa},par_T{ind_sin_sa}),2*pi);  
            if sign(par_T{par_option_T_season(i)}(1))~=sign(A_a)
                A_a=-A_a;
                phi_a=wrapToPi(phi_a-pi);
            end   
            if sign(par_T{par_option_T_season(i)}(1))~=sign(A_sa)
                A_sa=-A_sa;
                phi_sa=wrapToPi(phi_sa-pi);
            end               
            par_T{par_option_T_season(i)}(2:5)=[A_a,phi_a,A_sa,phi_sa];
            ind_base=ind_cos_sa;
        end
    end
end

if ~isempty(par_option_h_season)
    ind_base=length(par_option_h);
    for i=1:length(par_option_h_season)
        if par_option_h(par_option_h_season(i))==3 % Annual
            ind_sin=ind_base+1;
            ind_cos=ind_base+2;
            A=sqrt(par_h{ind_sin}^2+par_h{ind_cos}^2);
            phi=mod(atan2(par_h{ind_cos},par_h{ind_sin}),2*pi);
            if sign(par_h{par_option_h_season(i)}(1))~=sign(A)
                A=-A;
                phi=wrapToPi(phi-pi);
            end             
            par_h{par_option_h_season(i)}(2:3)=[A,phi];
            ind_base=ind_cos;            
        elseif par_option_h(par_option_h_season(i))==5 % Annual + Semi-Annual
            ind_sin_a=ind_base+1;
            ind_cos_a=ind_base+2;
            ind_sin_sa=ind_base+3;
            ind_cos_sa=ind_base+4;         
            A_a=sqrt(par_h{ind_sin_a}^2+par_h{ind_cos_a}^2);
            phi_a=mod(atan2(par_h{ind_cos_a},par_h{ind_sin_a}),2*pi);
            A_sa=sqrt(par_h{ind_sin_sa}^2+par_h{ind_cos_sa}^2);
            phi_sa=mod(atan2(par_h{ind_cos_sa},par_h{ind_sin_sa}),2*pi);
            if sign(par_T{par_option_h_season(i)}(1))~=sign(A_a)
                A_a=-A_a;
                phi_a=wrapToPi(phi_a-pi);
            end   
            if sign(par_T{par_option_h_season(i)}(1))~=sign(A_sa)
                A_sa=-A_sa;
                phi_sa=wrapToPi(phi_sa-pi);
            end                  
            par_h{par_option_h_season(i)}(2:5)=[A_a,phi_a,A_sa,phi_sa];
            ind_base=ind_cos_sa;            
        end
    end
end

%% Noise Parameters
% Perform Fitting for Noise Part using the Linear Regression Residual 
res_T_norm=res_T*sqrt(dt);
res_h_norm=res_h*sqrt(dt);

if n_T==1&n_h==1&n_g==2 % White-Additive 
sigma_T=std(res_T_norm);
sigma_h=std(res_h_norm);

B=NaN;
m_T=NaN;
m_h=NaN;

par_noise=[sigma_T;sigma_h;B;m_T;m_h];

elseif n_T==0&n_h==0&n_g==2 % Red-Additive 
    
if fitting_option_red=="LR"    
res_Tdiff=diff(res_T)/dt; % do not need to use normalized res_T 
res_hdiff=diff(res_h)/dt;

res_Ts=res_T(1:end-1);
res_hs=res_h(1:end-1);

[par_T_sigma,~,res_T_sigma]=regress_std(double(res_Tdiff),double(res_Ts));
[par_h_sigma,~,res_h_sigma]=regress_std(double(res_hdiff),double(res_hs));

m_T=-par_T_sigma;
m_h=-par_h_sigma;

elseif fitting_option_red=="AR1" 
m_T=-log(corr(res_T(1:end-1),res_T(2:end)))/dt;
m_h=-log(corr(res_h(1:end-1),res_h(2:end)))/dt;

elseif fitting_option_red=="ARn"
% m_T part 
[r,lags] = xcorr(res_T,30/dt,'normalized');
ind=find(lags<0);
if ~isempty(ind)
r(ind)=[];
lags(ind)=[];
end 
ind=min(find(r<0))-1;
if ~isempty(ind)
r=r(1:ind);
lags=lags(1:ind);
end 
y=-log(r);
X=lags*dt;
X=transpose(X);
[m_T]=regress_std(y,X);

% m_h part 
[r,lags] = xcorr(res_h,30,'normalized');
ind=find(lags<0);
if ~isempty(ind)
r(ind)=[];
lags(ind)=[];
end
ind=min(find(r<0))-1;
if ~isempty(ind)
r=r(1:ind);
lags=lags(1:ind);
end 
y=-log(r);
X=lags*dt;
X=transpose(X);
[m_h]=regress_std(y,X);
end 

sigma_T=std(res_T);
sigma_h=std(res_h);

B=NaN;

par_noise=[sigma_T;sigma_h;B;m_T;m_h];

elseif n_T==1&n_h==1&(n_g==0|n_g==1)% White-Multiplicative
n=10*12;
for i=1:1e5 
    ind=randi([1,length(res_T_norm)],n,1);
    res_T_var(i,1)=var(res_T_norm(ind));
    if n_g==0
        T_var(i,1)=var(Ts(ind));
    elseif n_g==1
        T_var(i,1)=var(heaviside(Ts(ind)).*Ts(ind));
    end 
end

par_T_res=regress_std(res_T_var,[T_var ones(length(T_var),1)]);

if par_T_res(2)>0&par_T_res(1)>0
    sigma_T=real(sqrt(par_T_res(2)));
    B=real(sqrt(par_T_res(1)/par_T_res(2))); 
else
    sigma_T=std(res_T_norm);
    B=0;     
end 

sigma_h=std(res_h_norm);
m_T=NaN;
m_h=NaN;

par_noise=[sigma_T;sigma_h;B;m_T;m_h];

elseif n_T==0&n_h==0&(n_g==0|n_g==1) % Red-Multiplicative
n=10*12;
for i=1:1e5 
    ind=randi([1,length(res_T_norm)],n,1);
    res_T_var(i,1)=var(res_T_norm(ind));
    if n_g==0
        T_var(i,1)=var(Ts(ind));
    elseif n_g==1
        T_var(i,1)=var(heaviside(Ts(ind)).*Ts(ind));
    end 
end

par_T_res=regress_std(res_T_var,[T_var ones(length(T_var),1)]);

if par_T_res(2)>0&par_T_res(1)>0
    sigma_T=real(sqrt(par_T_res(2)));
    B=real(sqrt(par_T_res(1)/par_T_res(2))); 
else
    sigma_T=std(res_T);
    B=0;     
end 
    
xi=res_T./(sigma_T*(1+B*Ts));
m_T=-regress_std(double(diff(xi)),double(xi(1:end-1)));

res_hdiff=diff(res_h)/dt;
res_hs=res_h(1:end-1);
[par_h_sigma,~,~]=regress_std(double(res_hdiff),double(res_hs));

m_h=-par_h_sigma;
sigma_h=std(res_h);

par_noise=[sigma_T;sigma_h;B;m_T;m_h];

end

par_noise=num2cell(par_noise);

%% Combine Parameters
par=[par_T(1:2);par_h(1:2);par_T(3:5);par_h(3);par_noise;n_T;n_h;n_g];

%% Using MAC method for Multiplicative T Noise Paremater Estimation 
if (n_g==0|n_g==1)&fitting_option_B=="MAC" 
par=RO_fitting_MAC(T,h,par); % Calculate sigma_T and B using MAC, and overwrite on the original par 
end 

end 

