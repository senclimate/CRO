function [T_out,h_out,noise_out]=RO_solver(par,IC,N,NE,NM,dt,saveat,savemethod,EF,noise_custom)

%% Setting Default Option
if nargin<5||isempty(NM)
    NM="EH";  
end
if nargin<6||isempty(dt)
    dt=0.1; 
end
if nargin<7||isempty(saveat)
    saveat=1.0;  
end
if nargin<8||isempty(savemethod)
    savemethod="sampling";
end
if nargin<9||isempty(EF)
    EF={0.0;0.0};  
end
if nargin<10||isempty(noise_custom)
    noise_custom=[];  
end

%% Setting Total Array Number
NT=round(N/dt);

%% Setting Output Variable 
T_out=zeros(NT,NE);
h_out=zeros(NT,NE);
noise_out=NaN(NT-1,4,NE);

%% Convert Ito to Stratonovich: Correction on R (only when white-multiplicaitve noise)
if NM=="EM" % apply when EM (based on Ito) apply Stratonovich -> Ito correction 
    if ~isempty(par{9}) & par{14}==1  % apply when B is given and white noise 
        if par{16}==0 % n_g==0: B*T linear multiplicative function
            par{1}(1)=par{1}(1)+0.5*((par{9}*par{11})^2); % R'=R+0.5*(sigma_T*B)^2
        elseif par{16}==1 % n_g==1: B*H(T)*T Heaviside multiplicative function
            par{1}(1)=par{1}(1)+0.25*((par{9}*par{11})^2); % R'=R+0.25*(sigma_T*B)^2
        end
    end
end

%% Setting Parameter & External Forcing
par_in=par_processing(par,N,dt,NT); % parameter input (total: 16)
EF_in=par_processing(EF,N,dt,NT); % external forcing input (total: 2)

if isempty(noise_custom)
    noise_custom=NaN(NT-1,4);
else
    if isscalar(noise_custom)
        rng(noise_custom)
        noise_custom=randn(NT-1,4);
    end 
end 

%% Solve RO
parfor i=1:NE
[T,h,noise]=RO_integral(par_in,EF_in,NM,NT,dt,IC,noise_custom); 
T_out(:,i)=T;
h_out(:,i)=h;
noise_out(:,:,i)=noise;
end 

%% Select Time Step % Save Noise Output
% note: noise_out does not follow this savemethod option 
if savemethod=="sampling"
    T_out=T_out(1:round(saveat/dt):end,:);
    h_out=h_out(1:round(saveat/dt):end,:);
elseif savemethod=="mean"
    T_out_sel=T_out(0.5*round(1/dt):end-0.5*round(1/dt)-1,:);
    h_out_sel=h_out(0.5*round(1/dt):end-0.5*round(1/dt)-1,:);    
    T_out_sel=reshape(T_out_sel,round(saveat/dt),[],size(T_out_sel,2));
    T_out_sel=squeeze(mean(T_out_sel,1));
    h_out_sel=reshape(h_out_sel,round(saveat/dt),[],size(h_out_sel,2));
    h_out_sel=squeeze(mean(h_out_sel,1));
    if NE>1
    T_out=[T_out(1,:);T_out_sel];
    h_out=[h_out(1,:);h_out_sel];
    elseif NE==1
    T_out=[T_out(1),T_out_sel];
    h_out=[h_out(1),h_out_sel];
    T_out=transpose(T_out);
    h_out=transpose(h_out);    
    end 
end

end 