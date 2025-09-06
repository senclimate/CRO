function par_out=par_processing(par_in,N,dt,NT)
% setting output
par_out=ones(NT,length(par_in));

% setting time
time=0:dt:(NT-1)*dt;  % time resolution of dt (month), length is NT
time_N=0:1:(N-1)+12; % time resolution of 1 (month), length is N

for i=1:length(par_in)
    
    % select par 
    par=par_in{i};
    
    % convert empty to zero
    if isempty(par)
        par=0.0;
    end 
    
    % convert vector to array (if input 'par' is vector)
    if size(par,1)>1
        par=transpose(par);
    end
    
    % convert par to format of [NTx1];
    if length(par)==1 % constant
        par_out(:,i)=par*ones(NT,1);
    elseif length(par)==NT % fine resolution input
        par_out(:,i)=par;
    elseif length(par)==N % intermediate resolution input
        par_interp=interp1(time_N,par,time);
        par_out(:,i)=par_interp;
    elseif length(par)==12 % seasonal value
        par_mat=repmat(par,[1,(N+12)/12]); % monthly
        par_interp=interp1(time_N,par_mat,time);
        par_out(:,i)=par_interp;
    elseif length(par)==5 % pre-defined sinusoidal cycle (annual + semi-annual)
        wa=(2*pi)/12;
        par_mat=par(1)+par(2)*sin(wa*time+par(3))+par(4)*sin(2*wa*time+par(5));
        par_out(:,i)=par_mat;
    elseif length(par)==3 % pre-defined sinusoidal cycle (annual)
        wa=(2*pi)/12;
        par_mat=par(1)+par(2)*sin(wa*time+par(3));
        par_out(:,i)=par_mat;      
    else
        display("Error: Please specifiy valid parameter value")
    end
    
end

end