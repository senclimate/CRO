function par=RO_fitting_MLE_white(T,h,par_option_T,par_option_h,par_option_noise,dt,varargin)
% Convert to Vector 
T=reshape(T,[length(T),1]);
h=reshape(h,[length(h),1]);

par_option_T=reshape(par_option_T,[1,length(par_option_T)]);
par_option_h=reshape(par_option_h,[1,length(par_option_h)]);
par_option_noise=reshape(par_option_noise,[1,length(par_option_noise)]);

n_T=par_option_noise(1);
n_h=par_option_noise(2);
n_g=par_option_noise(3);

n=length(T);

% Setting Option
par_option_T_onoff=par_option_T;
par_option_h_onoff=par_option_h;
par_option_T_onoff(find(par_option_T_onoff>=1))=1;
par_option_h_onoff(find(par_option_h_onoff>=1))=1;
par_option_T_h_onoff=[par_option_T_onoff par_option_h_onoff];

par_option_T_season=find(par_option_T>1);
par_option_h_season=find(par_option_h>1);
par_option_T_h_season=find([par_option_T par_option_h]>1);

par_option_T_h=[par_option_T par_option_h];

% Setting M Matrix 
MT=[T h T.^2 -T.^3 T.*h zeros(n,3)]*dt.*[par_option_T_onoff zeros(1,3)]; % dT/dt part 
Mh=[zeros(n,5) -T -h -T.^2]*dt.*[zeros(1,5) par_option_h_onoff]; % dh/dt part 

w=(2*pi)/12;
%t=1:dt:(n-1)*dt+1;
t=0:dt:(n-1)*dt;
t=transpose(t);

t_shift=t+0.5*dt; % for seasonal parameter estimation

if ~isempty(par_option_T_season)
    for i=1:length(par_option_T_season)
        if par_option_T(par_option_T_season(i))==3 % Annual
        MT_add=[MT(:,par_option_T_season(i)).*sin(w*t_shift),MT(:,par_option_T_season(i)).*cos(w*t_shift)];
        elseif par_option_T(par_option_T_season(i))==5 % Annual + Semi-Annual
        MT_add=[MT(:,par_option_T_season(i)).*sin(w*t_shift),MT(:,par_option_T_season(i)).*cos(w*t_shift),MT(:,par_option_T_season(i)).*sin(2*w*t_shift),MT(:,par_option_T_season(i)).*cos(2*w*t_shift)];      
        end 
        MT=[MT MT_add];
        Mh_add=zeros(size(MT_add));
        Mh=[Mh Mh_add];      
    end
end

if ~isempty(par_option_h_season)
    for i=1:length(par_option_h_season)
        if par_option_h(par_option_h_season(i))==3
        Mh_add=[Mh(:,5+par_option_h_season(i)).*sin(w*t_shift),Mh(:,5+par_option_h_season(i)).*cos(w*t_shift)];
        elseif par_option_h(par_option_h_season(i))==5
        Mh_add=[Mh(:,5+par_option_h_season(i)).*sin(w*t_shift),Mh(:,5+par_option_h_season(i)).*cos(w*t_shift),Mh(:,5+par_option_h_season(i)).*sin(2*w*t_shift),Mh(:,5+par_option_h_season(i)).*cos(2*w*t_shift)];    
        end 
        Mh=[Mh Mh_add];
        MT_add=zeros(size(Mh_add));
        MT=[MT MT_add];
    end
end

M=cat(3,MT,Mh);
M=permute(M,[3 2 1]);
Mtr=permute(M,[2 1 3]);

M=M(:,:,1:end-1);
Mtr=Mtr(:,:,1:end-1);

ind=find([par_option_T_onoff zeros(1,3)]==0&[zeros(1,5) par_option_h_onoff]==0);
M(:,ind,:)=[];
Mtr(ind,:,:)=[];

% Setting x diff 
x=transpose([diff(T) diff(h)]); 
x=permute(x,[1 3 2]); % making into 3D matrix 

% Initial Guess 
sigma_T=1.0; % Null Value - Not Practically Used
sigma_h=1.0; % Null Value - Not Practically Used
sigma_B=0.0; % Setting Zero 

% Setting Multiplicative Noise Factor
if n_g==0
    noise_g=ones(size(T)); % Linear (B*T)
elseif n_g==1 
    noise_g=heaviside(T); % Heavisdie (B*H(T)*T)
elseif n_g==2
    noise_g=0;
end

% Setting Iteration 
if n_g==0|n_g==1 % multiplicative 
    niter=100;
elseif n_g==2 % additive 
    niter=1;
end

for iter=1:niter
% Setting R Matrix 
RT=[(sigma_T+sigma_B*noise_g.*T).^2 zeros(n,1)]*dt; 
Rh=[zeros(n,1) ones(n,1)*sigma_h.^2]*dt;

R=cat(3,RT,Rh);
R=permute(R,[3 2 1]);

Rinv=1/R; % diagonal matrix 
Rinv(find(isinf(Rinv)))=0; 

R=R(:,:,1:end-1);
Rinv=Rinv(:,:,1:end-1);

% Calculate thetao
thetao_A=sum(pagemtimes(pagemtimes(Mtr,Rinv),M),3);
thetao_B=sum(pagemtimes(pagemtimes(Mtr,Rinv),x),3);
thetao=(thetao_A)^(-1)*thetao_B;

% Calculate R Matrix
R_A=x-pagemtimes(M,thetao);
R_B=permute(R_A,[2 1 3]);

% Calculate Noise for Additive Noise RO
R=sum(pagemtimes(R_A,R_B),3)/(n-1);

% Calculate Noise for Multiplicative Noise RO
if iter==1
    sigma_T=(R(1,1)/dt)^0.5; % Initial Guess for sigma_T using Additive Noise
    sigma_B=0.0; % Initial Guess for sigma_B
end

func=@(sigma_x)double([sum((squeeze(pagemtimes(R_A(1,:,:),R_B(:,1,:)))-dt*(sigma_x(1)+sigma_x(2)*noise_g(1:end-1).*T(1:end-1)).^2).*((sigma_x(1)+sigma_x(2)*noise_g(1:end-1).*T(1:end-1)).^(-3)).*T(1:end-1))/n;...
    sum((squeeze(pagemtimes(R_A(1,:,:),R_B(:,1,:)))-dt*(sigma_x(1)+sigma_x(2)*noise_g(1:end-1).*T(1:end-1)).^2).*((sigma_x(1)+sigma_x(2)*noise_g(1:end-1).*T(1:end-1)).^(-3)))/n]);

if verLessThan('matlab', '8.3')  
    options = optimset('Display', 'off');  % Use optimset instead of optimoptions
else
    options = optimoptions('fsolve', 'Display', 'off');
end
sigma_x_out=fsolve(func,[double(sigma_T);double(sigma_B)],options);
sigma_T=sigma_x_out(1); 
sigma_B=sigma_x_out(2);
sigma_h=(R(2,2)/dt)^0.5; 

end 

% Post-Processing Output
thetao_ann=thetao(1:length(find(par_option_T_h_onoff==1)));
thetao_season=thetao(length(find(par_option_T_h_onoff==1))+1:end);

par_T_h=zeros(8,1);
par_T_h(find(par_option_T_h_onoff==1))=thetao_ann; 
par_T_h=num2cell(par_T_h);

if ~isempty(par_option_T_h_season)
    ind_base=0;
    for i=1:length(par_option_T_h_season)
        if par_option_T_h(par_option_T_h_season(i))==3
        ind_sin=ind_base+1;
        ind_cos=ind_base+2;
        A=sqrt(thetao_season(ind_sin)^2+thetao_season(ind_cos)^2);
        phi=mod(atan2(thetao_season(ind_cos),thetao_season(ind_sin)),2*pi);
        if sign(par_T_h{par_option_T_h_season(i)}(1))~=sign(A)
            A=-A;
            phi=wrapToPi(phi-pi);
        end
        if isempty(varargin)
            par_T_h{par_option_T_h_season(i)}(2:3)=[A,phi];
        else
            par_T_h{par_option_T_h_season(i)}(2:3)=[thetao_season(ind_sin),thetao_season(ind_cos)];
        end
        ind_base=ind_cos;
        elseif par_option_T_h(par_option_T_h_season(i))==5
        ind_sin_a=ind_base+1;
        ind_cos_a=ind_base+2;
        ind_sin_sa=ind_base+3;
        ind_cos_sa=ind_base+4;        
        A_a=sqrt(thetao_season(ind_sin_a)^2+thetao_season(ind_cos_a)^2);
        phi_a=mod(atan2(thetao_season(ind_cos_a),thetao_season(ind_sin_a)),2*pi);
        A_sa=sqrt(thetao_season(ind_sin_sa)^2+thetao_season(ind_cos_sa)^2);
        phi_sa=mod(atan2(thetao_season(ind_cos_sa),thetao_season(ind_sin_sa)),2*pi);
        if sign(par_T_h{par_option_T_h_season(i)}(1))~=sign(A_a)
             A_a=-A_a;
             phi_a=wrapToPi(phi_a-pi);
        end   
        if sign(par_T_h{par_option_T_h_season(i)}(1))~=sign(A_sa)
             A_sa=-A_sa;
             phi_sa=wrapToPi(phi_sa-pi);
        end         
        if isempty(varargin)        
            par_T_h{par_option_T_h_season(i)}(2:5)=[A_a,phi_a,A_sa,phi_sa];
        else
            par_T_h{par_option_T_h_season(i)}(2:5)=[thetao_season(ind_sin_a),thetao_season(ind_cos_a),thetao_season(ind_sin_sa),thetao_season(ind_cos_sa)];
        end 
        ind_base=ind_cos_sa;        
        end
    end
end

par_T=par_T_h(1:5);
par_h=par_T_h(6:8); 

par_det=[par_T(1:2);par_h(1:2);par_T(3:5);par_h(3)];
par_noise=num2cell([sigma_T;sigma_h;sigma_B/sigma_T;NaN;NaN;n_T;n_h;n_g]);

par=[par_det;par_noise];