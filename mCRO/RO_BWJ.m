function BWJ=RO_BWJ(par)
%% Calculate BWJ index from Par 
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

%% Calculate BWJ
gr=(R-epsilon)/2;
w=sqrt(4*F1*F2-(R+epsilon)^2)/2;
BWJ=gr+j*w;

end 