function p=uhopftref(p,fac)
% uhopftref: uniform mesh refinement in t for hopf ; 
% p=uhopftref(p,fac)
% fac=increase/decrease factor for number of time-slices
ho=p.hopf; t=ho.t; tl1=round(fac*ho.tl); tln=tl1+5-rem(tl1,5); 
fprintf('old/new # of t grid points: %i, %i\n',ho.tl, tln); 
tn=linspace(0,1,tln); ho.t=tn; % new grid 
ynew=interp1(t,ho.y',ho.t); ho.y=ynew'; % interpolate relevant fields to new grid
ydnew=interp1(t,ho.y0d',ho.t); ho.y0d=ydnew'; 
tau=reshape(ho.tau(1:end-2-p.hopf.nqh),p.nu,ho.tl); % tangent 
tau=interp1(t,tau',ho.t); p.hopf.tl=tln; 
tau=reshape(tau',1,p.nu*tln); 
p.hopf.xi=p.hopf.xi/fac; % adapt weight of u in arclength
tau=[tau ho.tau(end-1-p.hopf.nqh:end)]; tau=tau/honorm(p,tau); 
ho.tl=tln; ho.tau=tau; p.hopf=ho; % update p.hopf 
end