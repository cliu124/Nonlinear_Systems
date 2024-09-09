function [L,Dphi,np,xx,yy,r]=CFLapNBC(p)
% CFlap: Laplacian for Chebychev-Fourier on (unit)disk
%
% also generates the discretization by p.nr points in r
%      and p.na points in angle and 
% returns (cartesian) grids xx,yy and r-point in r (to be used in dchebint) 
N=2*p.nr; M=p.na; g=[0 1 0; 0 1 0]; N2=N/2; 
[r,D2t,D1t]=cheb2bc(N,g); % phim=0, phip=0 for NBCs 
D1=D2t(1:N2,1:N2); D2=D2t(1:N2,N:-1:N2+1);
E1=D1t(1:N2,1:N2); E2=D1t(1:N2,N:-1:N2+1);
% t=theta coordinate, ranging from 0 to 2*pi (M must be even)
dt=2*pi/M; t=dt*(1:M)'; M2=M/2;
D2t=toeplitz([-pi^2/(3*dt^2)-1/6,.5*(-1).^(2:M)./sin(dt*(1:M-1)/2).^2]);
% Laplacian in polar coordinates:
R=diag(1./r(1:N2));  Z=zeros(M2); I=eye(M2);
L=kron(D1+R*E1,eye(M))+kron(D2+R*E2,[Z I;I Z])+kron(R^2,D2t);
R2=diag(r(1:N2)==1); % to set BCs, identify bdr points and set the rows 
L2=kron(R2*E1,eye(M)); idx=find(L2~=0); du=kron(R^2,D2t); % to D1t
% replace L at r=1 by D2t (\pa_r^2 with right BCs) 
L(idx)=du(idx); 
[rr,tt]=meshgrid(r(1:N2),[0; t]); [xx,yy]=pol2cart(tt,rr);
np=(size(xx,1)-1)*size(xx,2); 
column=[0 .5*(-1).^(1:M-1).*cot((1:M-1)*dt/2)]'; % for Dphi 
Dp=toeplitz(column,column([1 M:-1:2]));
Dphi=kron(diag(ones(1,N2)),Dp); 