function [f,fp]=fif(par,v)
bet=par(6); nlsw=par(7); 
if nlsw==1; f=bet*v.*(1-v); fp=bet*(1-2*v); % logistic
else f=-v.*(v-bet).*(v-1); fp=-bet+2*(bet+1)*v-3*v.^2; % bistab
end

