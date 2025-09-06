function par=RO_fitting(T,h,par_option_T,par_option_h,par_option_noise,method_fitting,dt)

%% Convert Parameter Options
par_option_T=cell2mat(struct2cell(par_option_T));
par_option_h=cell2mat(struct2cell(par_option_h));
par_option_noise=string(struct2cell(par_option_noise));

%% Setting Default Option
if nargin<6||isempty(method_fitting)
    method_fitting=func_default_fitting_method(par_option_T,par_option_h,par_option_noise);  
end

if nargin<7||isempty(dt)
    dt=1.0;
end

%% Convert Parameter Options
if par_option_noise(1)=="white"
    par_option_noise_array(1)=1;
elseif par_option_noise(1)=="red"
    par_option_noise_array(1)=0;
end

if par_option_noise(2)=="white"
    par_option_noise_array(2)=1;
elseif par_option_noise(2)=="red"
    par_option_noise_array(2)=0;
end

if par_option_noise(3)=="additive"
    par_option_noise_array(3)=2;
elseif par_option_noise(3)=="multiplicative"
    par_option_noise_array(3)=0;
elseif par_option_noise(3)=="multiplicative-H"
    par_option_noise_array(3)=1;
end

par_option_noise=par_option_noise_array;

%% Setting Implicit Fitting Method Option (internally specified in this function)
% estimation methods for red noise damping coefficient from linear regression residual
fitting_option_red="AR1"; % options: "LR" or "AR1" or "ARn" (defualt: "LR")

%% Perform Fitting 
if method_fitting=="LR-F" % Ito based 
    par=RO_fitting_LR(T,h,par_option_T,par_option_h,par_option_noise,dt,"F","LR",fitting_option_red);
elseif method_fitting=="LR-C" % Stratonovich based 
    par=RO_fitting_LR(T,h,par_option_T,par_option_h,par_option_noise,dt,"C","LR",fitting_option_red);
elseif method_fitting=="LR-F-MAC" % Ito based 
    par=RO_fitting_LR(T,h,par_option_T,par_option_h,par_option_noise,dt,"F","MAC",fitting_option_red);
elseif method_fitting=="LR-C-MAC" % Ito based 
    par=RO_fitting_LR(T,h,par_option_T,par_option_h,par_option_noise,dt,"C","MAC",fitting_option_red);
elseif method_fitting=="MLE" % Ito based 
    par=RO_fitting_MLE(T,h,par_option_T,par_option_h,par_option_noise,dt);   
end

%disp(sprintf("Fitting Method: %s",method_fitting));

%% Convert Ito to Stratonovich: Correction on R (only when white multiplicaitve noise)
if method_fitting~="LR-C" & par_option_T(1)~=0 & par_option_noise(3)~=2 & par_option_noise(1)==1 % apply correction (except LR-C) & par_option_T is turned on (1 3 5). & white noise
    if par{16}==0 % n_g==0: B*T linear multiplicative function 
        par{1}(1)=par{1}(1)-0.5*((par{9}*par{11})^2); % R'=R-0.5*(sigma_T*B)^2
    elseif par{16}==1 % n_g==1: B*H(T)*T Heaviside multiplicative function 
        par{1}(1)=par{1}(1)-0.25*((par{9}*par{11})^2); % R'=R-0.25*(sigma_T*B)^2
    end
end

%% Post Processing for Output
% Convert NaN to empty (NaN is not assgied parameter for fitting)
for k = 1:numel(par)
    x = par{k};
    if isnumeric(x)
        if isscalar(x) && (isnan(x) | x==0)
            par{k} = [];
        else
            par{k}(isnan(x)) = [];
        end
    end
end

end 


