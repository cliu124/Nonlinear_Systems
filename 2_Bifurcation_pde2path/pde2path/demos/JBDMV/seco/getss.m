function [as,us]=getss(par)
j0=par(1); T=par(3); as=j0/(T*(j0^2+1)); us=as+j0; 