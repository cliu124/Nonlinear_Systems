function xin=xinorm(z,xi,nq,xiq)
% XINORM: compute the weighted arclength-continuation xi-norm of z
%
%  xin=xinorm(z,xi,nq,xiq)
%  z=given vector, xi=value of xi, nq=number of aux. eqn., xiq=weight for aux. eqn.
len=length(z); ul2=norm(z(1:len-1-nq)); q=norm(z(len-nq:len)); lamd=z(len); 
xin=sqrt(xi*ul2^2+xiq*q^2+(1-(xi+xiq)/2)*lamd^2);
