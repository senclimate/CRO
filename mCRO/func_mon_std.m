function X_mon_std=func_mon_std(X,dt)
%dt
%if matches("monthly",dt)
    dt=1.0;
%end 

X=reshape(X,[length(X),1]);

X_mon=reshape(X,12/dt,size(X,1)/(12/dt)); % calendar month x year x ensemble 
X_mon_std=std(X_mon,[],2);

end 
