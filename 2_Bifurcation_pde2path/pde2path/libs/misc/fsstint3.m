function [p,t1,ts,nc]=fsstint3(p,t0,ts,dt,dtplot,nt,nc,N,del,varargin)
% fsstint: Fourier split step for NLS/GP 1D 
% INPUT: p: pde2path standard struct, with IC in p.u 
%           p.prep=1  for initial call which may mod p.u via del  
% t0,dt: start time, numerical time step 
% dtplot: time-steps for plots 
% nt: number of time steps 
% nc: global counter of outputs (for repeated calls) 
% N : number of Fourier modes
% del:  IC of type u=(1+del)*p.u 
% nf=varargin:  do boundary filtering of radiation each nf'th step
% 
% OUTPUT:  p
% t1: endtime
% ts: time-series of t,N,H 
% nc: global counter 
if nargin>9; nf=varargin{1}; else nf=0; end 
if nargin>10; xf=varargin{2}; else xf=1; end 
xg=getpte(p); x1=xf*xg(1); x2=xf*xg(end); L=x2-x1; dx=L/N; np=p.np; t=t0; 
x=linspace(-L/2,L/2-dx,N)'; k_vec=(2*pi/L)*[0:N/2 (-N/2+1):-1]'; % grids
try prep=p.prep; catch; prep=0; end 
if prep==1; % first call, starting from the REAL ground state
   u=interp1(xg,p.u(1:np),x,'linear','extrap'); p.prep=0; 
   save([p.file.dir '/0.mat'], 'x','u','t'); 
else  u=p.uc(1:N); 
end 
par=p.u(p.nu+1:end); lam=par(1); s=par(2); sig=par(3); ga=par(4); 
V=pot(x,s)+lam; phi=u; figure(3);  plot(x,real(u)); %,x,imag(u));  pause 
isteps=ceil(dtplot/dt); nsaves=floor(nt/isteps); shift=0; 
%define a filter for removing high k radiation 
sl=2; rad_filt=(tanh(sl*(x(x<0)-9*(-L/2+shift)/10))+1)/2;
rad_filt=[rad_filt; (-tanh(sl*(x(x>=0)-9*(L/2+shift)/10))+1)/2];
%figure(1); clf; plot(x,rad_filt); pause 
energy=dx*trapz(abs(u).^2); der_term=abs(ifft(1i*k_vec.*fft(u))).^2;
NL_term=ga*abs(u).^(2*sig+2)/(sig+1)+(V-lam).*abs(u).^2; 
h1=dx*trapz(der_term); h2=dx*trapz(NL_term); ham=h1+h2; 
ts=[ts, [t;energy;ham]]; 
for m=nc+1:nc+nsaves %TIME STEPPING      
    u=fft(u); u=exp(-1i*k_vec.^2*dt/2).*u; u=ifft(u); % iu_t + u_xx =0    
    for j=1:isteps-1                
        u=u.*exp(-1i*(ga*abs(u).^(2*sig)+V)*dt);    %the part: iu_t=Vu+ga|u|^(2sig)u =0            
        u=fft(u); u=exp(-1i*k_vec.^2*dt).*u; u=ifft(u); % iu_t + u_xx =0     
        if nf>0; if(mod(j,nf)==0); u=u.*rad_filt; end; end  % filtering 
        t=t+dt;
    end
    u=u.*exp(-1i*(ga*abs(u).^(2*sig)+V)*dt);    %the part: iu_t=Vu+ga|u|^(2sig)u =0        
    u=fft(u); u=exp(-1i*k_vec.^2*dt/2).*u; u=ifft(u); t=t+dt;   % iu_t + u_xx =0                   
    %conserved quantities
    energy=dx*trapz(abs(u).^2);  der_term=abs(ifft(1i*k_vec.*fft(u))).^2;
    NL_term=ga*abs(u).^(2*sig+2)/(sig+1)+(V-lam).*abs(u).^2; 
    ham=dx*trapz(der_term+NL_term);    ts=[ts, [t;energy;ham]]; 
    figure(1); plot(x,abs(u),'linewidth',2); hold on;
    plot(x,abs(phi),'--','linewidth',2); pause(0.1); hold off;     
    fname=[p.file.dir '/' sprintf('%i',m),'.mat']; save(fname,'u','t','ts'); 
end
t1=t; nc=m; p.uc=u; % store last soln as new IC 