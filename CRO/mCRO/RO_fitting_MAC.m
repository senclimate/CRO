function par_out=RO_fitting_MAC(T,h,par)
% post-processing par
for i=1:length(par)
    par_sel(i)=par{i}(1); % only select annual mean value
end 

% assign deterministic parameters 
R=par_sel(1); % this R is Ito sense directly getting from the linear regression
F1=par_sel(2);
F2=par_sel(3);
epsilon=par_sel(4);
b_T=par_sel(5);
c_T=par_sel(6);
d_T=par_sel(7);
b_h=par_sel(8);
n_g=par_sel(16);

% remove mean
T=T-mean(T);
h=h-mean(h);

if n_g==0 % linear function 
T2=mean(T.^2);
T3=mean(T.^3);
T4=mean(T.^4);
T5=mean(T.^5);
T6=mean(T.^6);
hT=mean(h.*T);
hT2=mean(h.*(T.^2));
hT3=mean(h.*(T.^3));
hT4=mean(h.*(T.^4));

TA=R*T2+F1*hT+b_T*T3-c_T*T4+d_T*hT2;
T2A=R*T3+F1*hT2+b_T*T4-b_T*(T2)^2-c_T*T5+c_T*T3*T2+d_T*hT3-d_T*hT*T2;
T3A=R*T4+F1*hT3+b_T*T5-b_T*T3*T2-c_T*T6+c_T*(T3^2)+d_T*hT4-d_T*hT*T3;

k=-(2/3)*(T3A-1.5*T2A*T3/T2-3*TA*T2)/(T4-(T3^2)/T2-(T2)^2); % k=(sigma_T*B)^2

sigma_T=(-2*TA-k*T2)^0.5;
B=-(T2A+k*T3)/(2*(sigma_T^2)*T2);

par_out=par;
par_out{9}=sigma_T;
par_out{11}=B;

elseif n_g==1 % Heaviside-linear function
T2=mean(T.^2);
T3=mean(T.^3);
T4=mean(T.^4);
T5=mean(T.^5);
T6=mean(T.^6);
hT=mean(h.*T);
hT2=mean(h.*(T.^2));
hT3=mean(h.*(T.^3));
hT4=mean(h.*(T.^4));

TA=R*T2+F1*hT+b_T*T3-c_T*T4+d_T*hT2;
T2A=R*T3+F1*hT2+b_T*T4-b_T*(T2)^2-c_T*T5+c_T*T3*T2+d_T*hT3-d_T*hT*T2;
T3A=R*T4+F1*hT3+b_T*T5-b_T*T3*T2-c_T*T6+c_T*(T3^2)+d_T*hT4-d_T*hT*T3;

T1p=T;
T2p=T.^2;
T3p=T.^3;
T4p=T.^4; 

T1p(find(T1p<0))=[];
T2p(find(T2p<0))=[];
T3p(find(T3p<0))=[];
T4p(find(T4p<0))=[];

T1p=mean(T1p);
T2p=mean(T2p);
T3p=mean(T3p);
T4p=mean(T4p);

k1=(4*TA*T2+T2A*(2*T3p/T2p-T1p*T2/T2p)-(4/3)*T3A)/(2*T4p-2*T2p*T2-2*(T3p^2)/T2p+T2*T3p*T1p/T2p); % k1=(sigma_T*B)^2
k2=-(1/(2*T2p))*(T2A+T3p*k1);  % k2=(sigma_T^2)*B

B=k1/k2;
sigma_T=(k2/B)^0.5;

par_out=par;
par_out{9}=sigma_T;
par_out{11}=B;

end 