function [xs,noise]=EM_scheme(x,f,g,dt,noise)
    if isnan(noise)
        noise=randn;
    end 
    dW=sqrt(dt)*noise;
    xs=x+f*dt+g*dW; 
end 