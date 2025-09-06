function [T,h,noise_out]=RO_integral(par,EF,NM,NT,dt,IC,noise_in); 

%% Setting Varible Array
T=ones(NT,1); % T 
h=ones(NT,1); % h 
xi_T=ones(NT,1); % xi_T 
xi_h=ones(NT,1); % xi_h

T(1)=IC(1);
h(1)=IC(2); 
xi_T(1)=0.0;
xi_h(1)=0.0;

noise_out=NaN(NT-1,4);

%% Solve RO 
if NM=="EM" % Euler-Maruyama Method  
    for i=1:NT-1
        [f_T,g_T,f_h,g_h,f_xi_T,g_xi_T,f_xi_h,g_xi_h]=RO_tendency(par(i,:),T(i),h(i),xi_T(i),xi_h(i),EF(i,:));
        [T(i+1),noise_out(i,1)]=EM_scheme(T(i),f_T,g_T,dt,noise_in(i,1));
        [h(i+1),noise_out(i,2)]=EM_scheme(h(i),f_h,g_h,dt,noise_in(i,2));
        [xi_T(i+1),noise_out(i,3)]=EM_scheme(xi_T(i),f_xi_T,g_xi_T,dt,noise_in(i,3));
        [xi_h(i+1),noise_out(i,4)]=EM_scheme(xi_h(i),f_xi_h,g_xi_h,dt,noise_in(i,4));
    end
end

if NM=="EH" % Euler-Heun Method (default)
    for i=1:NT-1
        [f_T,g_T,f_h,g_h,f_xi_T,g_xi_T,f_xi_h,g_xi_h]=RO_tendency(par(i,:),T(i),h(i),xi_T(i),xi_h(i),EF(i,:));
        % calculate intermediate value        
        [T_s,noise_out(i,1)]=EM_scheme(T(i),f_T,g_T,dt,noise_in(i,1));
        [h_s,noise_out(i,2)]=EM_scheme(h(i),f_h,g_h,dt,noise_in(i,2));
        [xi_T_s,noise_out(i,3)]=EM_scheme(xi_T(i),f_xi_T,g_xi_T,dt,noise_in(i,3));
        [xi_h_s,noise_out(i,4)]=EM_scheme(xi_h(i),f_xi_h,g_xi_h,dt,noise_in(i,4));
        [f_T_s,g_T_s,f_h_s,g_h_s,f_xi_T_s,g_xi_T_s,f_xi_h_s,g_xi_h_s]=RO_tendency(par(i+1,:),T_s,h_s,xi_T_s,xi_h_s,EF(i+1,:));
        % calculate next step value           
        [T(i+1),~]=EM_scheme(T(i),0.5*(f_T+f_T_s),0.5*(g_T+g_T_s),dt,noise_out(i,1));
        [h(i+1),~]=EM_scheme(h(i),0.5*(f_h+f_h_s),0.5*(g_h+g_h_s),dt,noise_out(i,2));
        [xi_T(i+1),~]=EM_scheme(xi_T(i),0.5*(f_xi_T+f_xi_T_s),0.5*(g_xi_T+g_xi_T_s),dt,noise_out(i,3));
        [xi_h(i+1),~]=EM_scheme(xi_h(i),0.5*(f_xi_h+f_xi_h_s),0.5*(g_xi_h+g_xi_h_s),dt,noise_out(i,4));     
    end
end

end