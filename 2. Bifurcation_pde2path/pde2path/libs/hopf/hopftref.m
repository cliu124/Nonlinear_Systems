function p=hopftref(p,t0)
% hopftref: mesh refinement in t for hopf; 
% p=hopftref(p,t0)
% time-slice(s) t0 BEHIND which to refine 
% (bisection of the next 5 intervals) given by user; 
ho=p.hopf; t=ho.t;
[tv,idx]=sort(abs(p.hopf.T*t(:)-t0)); ti=idx(1);  ti, p.hopf.T*t(ti+2)
if ti>ho.tl-5; ti=ho.tl-5; end 
t1=t(1:ti); t2=t(ti+1:end); 
t2t=[(t1(end)+t2(1))/2 t2(1) (t2(1)+t2(2))/2 t2(2) (t2(2)+t2(3))/2 ... 
    t2(3) (t2(3)+t2(4))/2 t2(4) (t2(4)+t2(5))/2 t2(5:end)]; 
ho.t=[t1 t2t]; % new grid 
ynew=interp1(t,ho.y',ho.t); ho.y=ynew'; % interpolate relevant fields to new grid
if 1; ydnew=interp1(t,ho.y0d',ho.t); ho.y0d=ydnew'; 
else par=p.u(p.nu+1:end);
for i=1:ho.tl+5; % phase cond 
      f=hosrhs(0,ho.y(:,i),p,par,ho.T); 
      ho.y0d(1:p.nu,i)=f(1:p.nu); 
end; 
end 
tau=reshape(ho.tau(1:end-2-p.hopf.nqh),p.nu,ho.tl); % tangent 
tau=interp1(t,tau',ho.t); 
tau=reshape(tau',1,p.nu*(ho.tl+5)); 
tau=[tau ho.tau(end-1-p.hopf.nqh:end)]; tau=tau/honorm(p,tau); 
ho.tau=tau; ho.tl=ho.tl+5; p.hopf=ho; % update p.hopf 
end 
