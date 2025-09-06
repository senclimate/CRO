function fitting_method=func_default_fitting_method(par_option_T,par_option_h,par_option_noise)

%% Setting Method 
if any(par_option_T>1)|any(par_option_h>1)
    seasonal_type="seasonal";
else
    seasonal_type="constant";
end 

if par_option_T(3)==0&par_option_T(4)==0&par_option_T(5)==0&par_option_h(3)==0
    det_type="linear";
else
    det_type="nonlinear";
end 

noise_color_type=par_option_noise(1);
noise_amp_type=par_option_noise(3);

%% Load Table and Find Method
method_table=readtable('table_default_fitting_method.xlsx');

is_match = strcmp(method_table.seasonal_type, seasonal_type) & ...
           strcmp(method_table.det_type, det_type) & ...
           strcmp(method_table.noise_color_type, noise_color_type) & ...
           strcmp(method_table.noise_amp_type, noise_amp_type);
       
fitting_method = string(method_table.fitting_method(is_match));       
end 

