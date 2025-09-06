function par=RO_fitting_MLE(T,h,par_option_T,par_option_h,par_option_noise,dt)

n_T=par_option_noise(1);
n_h=par_option_noise(2);

if n_T==1&n_h==1 % white noise RO
    par=RO_fitting_MLE_white(T,h,par_option_T,par_option_h,par_option_noise,dt);
elseif n_T==0&n_h==0 % red noise RO
    par=RO_fitting_MLE_red(T,h,par_option_T,par_option_h,par_option_noise,dt);
else
    disp("Error!")
end 