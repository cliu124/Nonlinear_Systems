function [L,Dphi,np,r,t]=CFLapDBC(p)
% CFlap: Laplacian for Chebychev-Fourier on (unit)disk, 
% to act on full mesh, where u in L*u must first be 'filled' by BCvalues 
%
% also generates the dicretization by p.nr points in r
%      and p.na points in angle in t 
% Following Trefethen, Ch.11 
N=2*p.nr-1; M=p.na; N2=(N-1)/2; [D,r]=cheb(N); D2=D^2; 
D1=D2(1:N2+1,1:N2+1); D2=D2(1:N2+1,N+1:-1:N2+2); 
E1=D(1:N2+1,1:N2+1); E2=D(1:N2+1,N+1:-1:N2+2);
% t=theta coordinate, ranging from 0 to 2*pi (M must be even)
dt=2*pi/M; t=dt*(1:M)'; M2=M/2;
D2t=toeplitz([-pi^2/(3*dt^2)-1/6,.5*(-1).^(2:M)./sin(dt*(1:M-1)/2).^2]);
% Laplacian in polar coordinates:
R=diag(1./r(1:N2+1));  Z=zeros(M2); I=eye(M2);
L=kron(D1+R*E1,eye(M))+kron(D2+R*E2,[Z I;I Z])+kron(R^2,D2t);
r=r(1:N2+1); np=(p.nr-1)*p.na;  
column=[0 .5*(-1).^(1:M-1).*cot((1:M-1)*dt/2)]'; % for Dphi 
Dp=toeplitz(column,column([1 M:-1:2]));
Dphi=kron(diag(ones(1,N2+1)),Dp); 