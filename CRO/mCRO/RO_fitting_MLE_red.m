function par=RO_fitting_mle_red(T,h,par_option_T,par_option_h,par_option_noise,dt)
% Convert to Vector 
T=reshape(T,[length(T),1]);
h=reshape(h,[length(h),1]);

par_option_T=reshape(par_option_T,[1,length(par_option_T)]);
par_option_h=reshape(par_option_h,[1,length(par_option_h)]);
par_option_noise=reshape(par_option_noise,[1,length(par_option_noise)]);

n_T=par_option_noise(1);
n_h=par_option_noise(2);
n_g=par_option_noise(3);

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

% Initial Guess for White Noise Fitting Value
par_white_in=RO_fitting_MLE_white(T,h,par_option_T,par_option_h,[1 1 par_option_noise(3)],dt,'raw');

for i=1:length(par_white_in)
    par_white(i)=par_white_in{i}(1); % only select annual mean value
end 

R=par_white(1);
F1=par_white(2);
F2=par_white(3);
epsilon=par_white(4);
b_T=par_white(5);
c_T=par_white(6);
d_T=par_white(7);
b_h=par_white(8);
sigma_T=par_white(9);
sigma_h=par_white(10);
sigma_B=par_white(11)*sigma_T;
m_T=1.0;
m_h=1.0;

% Setting Seasonal Parameters 
par_season=[];
par_season_T=[];
par_season_h=[];

par_order_array=[1 2 5 6 7 3 4 8]; % T (5) and h (3)

if ~isempty(par_option_T_h_season)    
    for i=1:length(par_order_array)
        if length(par_white_in{par_order_array(i)})>1
            par_season=[par_season par_white_in{par_order_array(i)}(2:end)];
            if i<=5
                par_season_T=[par_season_T par_white_in{par_order_array(i)}(2:end)];
            else
                par_season_h=[par_season_h par_white_in{par_order_array(i)}(2:end)];
            end
        end
    end
end

sigma_T2=0.01;
sigma_h2=0.01;

% Interpolation of T and h 
dt_interp=0.01; % empirically dt_interp should be 0.01 or smaller 
t=0:dt:(length(T)-1)*dt;

if dt>dt_interp
time_interp=0:dt_interp:(length(T)-1)*dt;

T_interp=interp1(transpose(t),T,transpose(time_interp));
h_interp=interp1(transpose(t),h,transpose(time_interp));

T=T_interp;
h=h_interp;
dt=dt_interp;
t=time_interp;
end 

w=(2*pi)/12;
t=transpose(t);

t_shift=t+0.5*dt; % for seasonal parameter estimation

n=length(T);

% Setting x diff 
x=transpose([diff(T) diff(h)]); 
x=permute(x,[1 3 2]); % making into 3D matrix 

% Setting noise factor
if n_g==0
    noise_g=ones(size(T)); % Linear (B*T)
elseif n_g==1 
    noise_g=heaviside(T); % Heavisdie (B*H(T)*T)
elseif n_g==2 % Additive (B=0)
    noise_g=zeros(size(T));
end

niter=100; % 200

for iter=1:niter
    
% Setting CGNS Matrix
A0_T=[R*T+F1*h+b_T*T.^2-c_T*T.^3+d_T*T.*h];
A0_h=[-F2*T-epsilon*h-b_h*T.^2];

Tv=[T,h,T.^2,-T.^3,T.*h];
hv=[-T,-h,-T.^2];

if ~isempty(par_option_T_season)
    ind_base=0;
    for i=1:length(par_option_T_season)
        if par_option_T(par_option_T_season(i))==3 % annual
            ind_sin=ind_base+1;
            ind_cos=ind_base+2;
            A0_T_add=par_season_T(ind_sin)*Tv(:,par_option_T_season(i)).*sin(w*t_shift)+par_season_T(ind_cos)*Tv(:,par_option_T_season(i)).*cos(w*t_shift);
            A0_T=A0_T+A0_T_add;
            ind_base=ind_cos; 
        elseif par_option_T(par_option_T_season(i))==5 % annual + semi-annual
            ind_sin_a=ind_base+1;
            ind_cos_a=ind_base+2;
            ind_sin_sa=ind_base+3;
            ind_cos_sa=ind_base+4;            
            A0_T_add_a=par_season_T(ind_sin_a)*Tv(:,par_option_T_season(i)).*sin(w*t_shift)+par_season_T(ind_cos_a)*Tv(:,par_option_T_season(i)).*cos(w*t_shift);
            A0_T_add_sa=par_season_T(ind_sin_sa)*Tv(:,par_option_T_season(i)).*sin(2*w*t_shift)+par_season_T(ind_cos_sa)*Tv(:,par_option_T_season(i)).*cos(2*w*t_shift);
            A0_T=A0_T+A0_T_add_a+A0_T_add_sa;
            ind_base=ind_cos_sa;             
        end 
    end
end

if ~isempty(par_option_h_season)
    ind_base=0;    
    for i=1:length(par_option_h_season)
        if par_option_h(par_option_h_season(i))==3 % annual
            ind_sin=ind_base+1;
            ind_cos=ind_base+2;
            A0_h_add=par_season_h(ind_sin)*Tv(:,par_option_h_season(i)).*sin(w*t_shift)+par_season_h(ind_cos)*Tv(:,par_option_h_season(i)).*cos(w*t_shift);
            A0_h=A0_h+A0_h_add;
            ind_base=ind_cos; 
        elseif par_option_h(par_option_h_season(i))==5 % annual + semi-annual
            ind_sin_a=ind_base+1;
            ind_cos_a=ind_base+2;
            ind_sin_sa=ind_base+3;
            ind_cos_sa=ind_base+4;            
            A0_h_add_a=par_season_h(ind_sin_a)*Tv(:,par_option_h_season(i)).*sin(w*t_shift)+par_season_h(ind_cos_a)*Tv(:,par_option_h_season(i)).*cos(w*t_shift);
            A0_h_add_sa=par_season_h(ind_sin_sa)*Tv(:,par_option_h_season(i)).*sin(2*w*t_shift)+par_season_h(ind_cos_sa)*Tv(:,par_option_h_season(i)).*cos(2*w*t_shift);
            A0_h=A0_h+A0_h_add_a+A0_h_add_sa;
            ind_base=ind_cos_sa;             
        end 
    end
end

A0=[A0_T A0_h];
A0=transpose(A0);

A1_T=[sigma_T+sigma_B*noise_g.*T zeros(size(T))];
A1_h=[zeros(size(T)) sigma_h*ones(size(T))];
A1=cat(3,A1_T,A1_h);
A1=permute(A1,[3 2 1]);

B1=[sigma_T2 0;0 sigma_h2];
a0=[0;0];
a1=[-m_T 0;0 -m_h];
b2=[1 0;0 1];

% Solve mu_f and R_f 
mu_f=zeros(2,length(T));
R_f=zeros(2,2,length(T));
R_f(:,:,1)=0.01*eye(2);

for i=1:n-1
A0_sel=squeeze(A0(:,i));
A1_sel=squeeze(A1(:,:,i));
R_f_sel=squeeze(R_f(:,:,i));
mu_f_sel=squeeze(mu_f(:,i));
x_sel=squeeze(x(:,:,i));

mu_f_tendency=a0+a1*mu_f(:,i)+(R_f_sel*transpose(A1_sel))*((B1*transpose(B1))^(-1))*(x_sel/dt-(A0_sel+A1_sel*mu_f_sel));
R_f_tendency=a1*R_f_sel+R_f_sel*transpose(a1)+b2*transpose(b2)-(R_f_sel*transpose(A1_sel))*((B1*transpose(B1))^(-1))*(A1_sel*R_f_sel);

mu_f(:,i+1)=mu_f(:,i)+mu_f_tendency*dt;
R_f(:,:,i+1)=R_f(:,:,i)+R_f_tendency*dt;
end 

% Solve mu_s and R_s
mu_s=zeros(2,length(T));
R_s=zeros(2,2,length(T));

mu_s(:,end)=mu_f(:,end);
R_s(:,:,end)=R_f(:,:,end);

for i=n-1:-1:1
R_f_sel=squeeze(R_f(:,:,i+1));
mu_f_sel=squeeze(mu_f(:,i+1));
R_s_sel=squeeze(R_s(:,:,i+1));
mu_s_sel=squeeze(mu_s(:,i+1));

mu_s_tendency=-a0-a1*mu_s_sel+(b2*transpose(b2))*(R_f_sel^(-1))*(mu_f_sel-mu_s_sel);
R_s_tendency=-(a1+(b2*transpose(b2))*(R_f_sel^(-1)))*R_s_sel-R_s_sel*(transpose(a1)+(b2*transpose(b2))*R_f_sel)+b2*transpose(b2);

mu_s(:,i)=mu_s(:,i+1)+mu_s_tendency*dt;
R_s(:,:,i)=R_s(:,:,i+1)+R_s_tendency*dt;
end 

% Calculate C matrix and y moments 
C=zeros(size(R_f));
for i=1:size(C,3)
R_f_sel=squeeze(R_f(:,:,i));
C(:,:,i)=R_f_sel*transpose(eye(2,2)+a1*dt)*((b2*transpose(b2)*dt+(eye(2,2)+a1*dt)*R_f_sel*transpose(eye(2,2)+a1*dt))^(-1));
end 

yi=mu_s;
yiyi=zeros(size(R_s));
for i=1:size(R_s,3)
    R_s_sel=squeeze(R_s(:,:,i));
    mu_s_sel=squeeze(mu_s(:,i));
    yiyi(:,:,i)=R_s_sel+mu_s_sel*transpose(mu_s_sel);
end 

yi1yi=zeros(size(R_s,1),size(R_s,2),size(R_s,3)-1);
for i=1:size(R_s,3)-1
    yi1yi(:,:,i)=squeeze(R_s(:,:,i+1))*transpose(squeeze(C(:,:,i)))+squeeze(mu_s(:,i+1))*transpose(squeeze(mu_s(:,i)));
end 

% Setting M and Rinv 
Rinv=diag([1/sigma_T2^2,1/sigma_h2^2,1,1])*(1/dt);
Rinv=repmat(Rinv,[1 1 n]);

xi_T=transpose(yi(1,:));
xi_h=transpose(yi(2,:));

MT=[T h T.^2 -T.^3 T.*h xi_T xi_T.*noise_g.*T zeros(n,6)]*dt.*[par_option_T_onoff ones(1,2) zeros(1,6)]; % dT/dt part 
Mh=[zeros(n,7) -T -h -T.^2 xi_h zeros(n,2)]*dt.*[zeros(1,7) par_option_h_onoff ones(1,1) zeros(1,2)]; % dh/dt part 
MxiT=[zeros(n,11) -xi_T zeros(n,1)]*dt.*[zeros(1,11) ones(1,1) zeros(1,1)]; % dh/dt part 
Mxih=[zeros(n,12) -xi_h]*dt.*[zeros(1,11) zeros(1,1) ones(1,1)]; % dh/dt part 

if ~isempty(par_option_T_season)
    for i=1:length(par_option_T_season)
        if par_option_T(par_option_T_season(i))==3
        MT_add=[MT(:,par_option_T_season(i)).*sin(w*t_shift),MT(:,par_option_T_season(i)).*cos(w*t_shift)];
        elseif par_option_T(par_option_T_season(i))==5
        MT_add=[MT(:,par_option_T_season(i)).*sin(w*t_shift),MT(:,par_option_T_season(i)).*cos(w*t_shift),MT(:,par_option_T_season(i)).*sin(2*w*t_shift),MT(:,par_option_T_season(i)).*cos(2*w*t_shift)];    
        end 
        MT=[MT MT_add];        
        Mh_add=zeros(size(MT_add));
        Mh=[Mh Mh_add];
        MxiT_add=zeros(size(MT_add));
        MxiT=[MxiT MxiT_add];
        Mxih_add=zeros(size(MT_add));
        Mxih=[Mxih Mxih_add];        
    end
end

if ~isempty(par_option_h_season)
    for i=1:length(par_option_h_season)
        if par_option_h(par_option_h_season(i))==3
        Mh_add=[Mh(:,7+par_option_h_season(i)).*sin(w*t_shift),Mh(:,7+par_option_h_season(i)).*cos(w*t_shift)];
        elseif par_option_h(par_option_h_season(i))==5
        Mh_add=[Mh(:,7+par_option_h_season(i)).*sin(w*t_shift),Mh(:,7+par_option_h_season(i)).*cos(w*t_shift),Mh(:,7+par_option_h_season(i)).*sin(2*w*t_shift),Mh(:,7+par_option_h_season(i)).*cos(2*w*t_shift)];    
        end 
        Mh=[Mh Mh_add];
        MT_add=zeros(size(Mh_add));
        MT=[MT MT_add];
        MxiT_add=zeros(size(Mh_add));
        MxiT=[MxiT MxiT_add];
        Mxih_add=zeros(size(Mh_add));
        Mxih=[Mxih Mxih_add];             
    end
end

M=cat(3,MT,Mh,MxiT,Mxih);
M=permute(M,[3 2 1]);
Mtr=permute(M,[2 1 3]);

if n_g==0|n_g==1 % multi
par_option_matrix=[par_option_T_onoff ones(1,1) ones(1,1) zeros(1,6);zeros(1,7) par_option_h_onoff ones(1,1) zeros(1,2);zeros(1,11) ones(1,1) zeros(1,1);zeros(1,11) zeros(1,1) ones(1,1)];
elseif n_g==2 % additive
par_option_matrix=[par_option_T_onoff ones(1,1) zeros(1,1) zeros(1,6);zeros(1,7) par_option_h_onoff ones(1,1) zeros(1,2);zeros(1,11) ones(1,1) zeros(1,1);zeros(1,11) zeros(1,1) ones(1,1)];    
end 

par_option_matrix=[par_option_matrix(1,:) ones(1,length(par_season_T)) zeros(1,length(par_season_h))...
    ;par_option_matrix(2,:) zeros(1,length(par_season_T)) ones(1,length(par_season_h))...
    ;par_option_matrix(3,:) zeros(1,length(par_season_T)+length(par_season_h))...
    ;par_option_matrix(4,:) zeros(1,length(par_season_T)+length(par_season_h))];

ind=find(sum(par_option_matrix,1)==0);
M(:,ind,:)=[];
Mtr(ind,:,:)=[];

% Calculate thetao_A
thetao_A=pagemtimes(pagemtimes(Mtr,Rinv),M);

nM=size(thetao_A,1)-length(par_season_T)-length(par_season_h);
nT=length(find(par_option_T_onoff==1));

thetao_A(nT+1,nT+1,:)=dt*squeeze(yiyi(1,1,:))/(sigma_T2.^2);
thetao_A(nM-2,nM-2,:)=dt*squeeze(yiyi(2,2,:))/(sigma_h2.^2);
thetao_A(nM-1,nM-1,:)=dt*squeeze(yiyi(1,1,:))/(1);
thetao_A(nM,nM,:)=dt*squeeze(yiyi(2,2,:))/(1);

if n_g==0|n_g==1
    thetao_A(nT+2,nT+2,:)=dt*(noise_g).*(T.^2).*squeeze(yiyi(1,1,:))/(sigma_T2.^2);
end

thetao_A=thetao_A(:,:,1:end-1);

% Calculate thetao_B
x_NaN=zeros(size(x));
x_all=cat(1,x,x_NaN);

thetao_B=pagemtimes(pagemtimes(Mtr(:,:,1:end-1),Rinv(:,:,1:end-1)),x_all);
thetao_B(nM-1,1,:)=(-squeeze(yi1yi(1,1,:))+squeeze(yiyi(1,1,1:end-1)))/(1);
thetao_B(nM,1,:)=(-squeeze(yi1yi(2,2,:))+squeeze(yiyi(2,2,1:end-1)))/(1);

% Calculate thetao
ind_res=1;

thetao_A_sum=sum(thetao_A(:,:,ind_res:end-ind_res+1),3);
thetao_B_sum=sum(thetao_B(:,:,ind_res:end-ind_res+1),3);
thetao=(thetao_A_sum)^(-1)*thetao_B_sum;

% Calculate R 
R_A=x_all-pagemtimes(M(:,:,1:end-1),thetao);
R_B=permute(R_A,[2 1 3]);

R_AB=pagemtimes(R_A,R_B);

R_AB(1,1,:)=squeeze(R_AB(1,1,:))-((sigma_T*xi_T(1:end-1,1)+sigma_B*noise_g(1:end-1,:).*xi_T(1:end-1,1).*T(1:end-1,1))*dt).^2;
R_AB(2,2,:)=squeeze(R_AB(2,2,:))-((sigma_h*xi_h(1:end-1,1))*dt).^2;

R_AB(1,1,:)=squeeze(R_AB(1,1,:))+(squeeze(yiyi(1,1,1:end-1))).*(dt*(sigma_T+sigma_B*noise_g(1:end-1).*T(1:end-1))).^2;
R_AB(2,2,:)=squeeze(R_AB(2,2,:))+(squeeze(yiyi(2,2,1:end-1))).*(dt*sigma_h).^2;

R_AB(3,3,:)=squeeze(yiyi(1,1,2:end))+squeeze(yiyi(1,1,1:end-1))*(thetao(end-1-length(par_season_T)-length(par_season_h))*dt-1)^2+2*squeeze(yi1yi(1,1,:))*(thetao(end-1-length(par_season_T)-length(par_season_h))*dt-1);
R_AB(4,4,:)=squeeze(yiyi(2,2,2:end))+squeeze(yiyi(2,2,1:end-1))*(thetao(end-length(par_season_T)-length(par_season_h))*dt-1)^2+2*squeeze(yi1yi(2,2,:))*(thetao(end-length(par_season_T)-length(par_season_h))*dt-1);

R_M=sum(R_AB(:,:,ind_res:end-ind_res+1),3)/(size(R_AB,3));

% Save to Parameter
thetao_out=zeros(13+length(par_season_T)+length(par_season_h),1);
thetao_out(find(sum(par_option_matrix,1)==1))=thetao; 
thetao_out_save(:,iter)=thetao_out;

% Convergence Acceleration
% if iter>=2&iter<=floor(niter)/2
% ap=1+1/iter;
% thetao_out=thetao_out_save(:,iter-1)+ap*(thetao_out_save(:,iter)-thetao_out_save(:,iter-1));
% end 

R=thetao_out(1);
F1=thetao_out(2);
b_T=thetao_out(3);
c_T=thetao_out(4);
d_T=thetao_out(5);
sigma_T=thetao_out(6);
sigma_B=thetao_out(7);
B=sigma_B/sigma_T;
F2=thetao_out(8);
epsilon=thetao_out(9);
b_h=thetao_out(10);
sigma_h=thetao_out(11);
m_T=thetao_out(12);
m_h=thetao_out(13);

if ~isempty(par_option_T_h_season)
    par_season=thetao_out(14:end);
    if ~isempty(par_option_T)
        par_season_T=thetao_out(14:13+length(par_season_T));
    end 
    if ~isempty(par_option_h)
        par_season_h=thetao_out(14+length(par_season_T):end);
    end     
end 

sigma_T2=0.01;
sigma_h2=0.01;

end 

% Setting Final Output
par=[R;F1;F2;epsilon;b_T;c_T;d_T;b_h;sigma_T/sqrt(2*m_T);sigma_h/sqrt(2*m_h);B;m_T;m_h;n_T;n_h;n_g];
par=num2cell(par);

% Setting Final Output
if ~isempty(par_option_T_h_season)
    ind_base=0;
    for i=1:length(par_option_T_h_season)
        if par_option_T_h(par_option_T_h_season(i))==3
        ind_sin=ind_base+1;
        ind_cos=ind_base+2;
        A=sqrt(par_season(ind_sin)^2+par_season(ind_cos)^2);
        phi=mod(atan2(par_season(ind_cos),par_season(ind_sin)),2*pi);
        if sign(par{par_option_T_h_season(i)}(1))~=sign(A)
            A=-A;
            phi=wrapToPi(phi-pi);
        end
        par{par_order_array(par_option_T_h_season(i))}(2:3)=[A,phi];
        ind_base=ind_cos;
        elseif par_option_T_h(par_option_T_h_season(i))==5
        ind_sin_a=ind_base+1;
        ind_cos_a=ind_base+2;
        ind_sin_sa=ind_base+3;
        ind_cos_sa=ind_base+4;        
        A_a=sqrt(par_season(ind_sin_a)^2+par_season(ind_cos_a)^2);
        phi_a=mod(atan2(par_season(ind_cos_a),par_season(ind_sin_a)),2*pi);
        A_sa=sqrt(par_season(ind_sin_sa)^2+par_season(ind_cos_sa)^2);
        phi_sa=mod(atan2(par_season(ind_cos_sa),par_season(ind_sin_sa)),2*pi);
        if sign(par{par_option_T_h_season(i)}(1))~=sign(A_a)
             A_a=-A_a;
             phi_a=wrapToPi(phi_a-pi);
        end   
        if sign(par{par_option_T_h_season(i)}(1))~=sign(A_sa)
             A_sa=-A_sa;
             phi_sa=wrapToPi(phi_sa-pi);
        end                
        par{par_order_array(par_option_T_h_season(i))}(2:5)=[A_a,phi_a,A_sa,phi_sa];
        ind_base=ind_cos_sa;        
        end
    end
end




