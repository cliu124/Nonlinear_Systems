function p=hoada(p,bis)
% hoada: mesh refinement in t for hopf; 
% p=hoada(p,bis)
% bis=number of intervals to bisect
ho=p.hopf; t=ho.t; tl=length(t); h=diff(t); y=ho.y; lam=y(p.nu,:); 
mclf(12); plot(t,lam,'-*'); xlabel('t'); title('\lambda_1'); grid on;  
try Xsw=p.sw.Xcont; catch; Xsw=0; end
Dt=makeDtm(t); Dt2=makeDtm2(t); l1=Dt*lam'; l2=Dt2*lam'; 
mclf(13);  plot(t,l2); xlabel('t'); title('\lambda_1´´'); grid on;  
[l2s,idx]=sort(abs(l2),'descend'); %%idx=setdiff(idx,[1 tl]); 
idx=idx(1:bis); [idx,ia]=setdiff(idx, [1 tl], 'stable'); 
bis=bis-(bis-length(ia)); 
[ids,id]=sort(idx,'ascend'); tn=t; 
for i=1:bis 
 ti=ids(i); t0=tn(ti); 
 t1=tn(1:ti-1); t2=tn(ti+1:end); % t1(end), t0, t2(1)
 t2n=[(t1(end)+t0)/2, t0, (t0+t2(1))/2]; 
 tn=[t1 t2n t2]; % new grid 
 ids(i+1:end)=ids(i+1:end)+2; % two additional points 
end 
tln=length(tn); 
ynew=interp1(ho.t,ho.y',tn,'spline'); ho.y=ynew'; % interpolate relevant fields to new grid
ydnew=interp1(ho.t,ho.y0d',tn); ho.y0d=ydnew'; 
tau=reshape(ho.tau(1:end-2-p.hopf.nqh),p.nu,ho.tl); % tangent 
tau=interp1(ho.t,tau',tn, 'spline'); 
tau=reshape(tau',1,p.nu*(ho.tl+2*bis)); 
tau=[tau ho.tau(end-1-p.hopf.nqh:end)]; tau=tau/honorm(p,tau); 
ho.tau=tau; ho.tl=length(tn); 
if Xsw>0
  X=ho.X; X1=X(:,1,:);  X1=reshape(X1,p.np,tl); 
  X2=X(:,2,:);  X2=reshape(X2,p.np,tl);
  X1n=interp1(ho.t,X1',tn,'spline'); X1n=X1n'; 
  X2n=interp1(ho.t,X2',tn, 'spline'); X2n=X2n'; 
  for i=1:tln 
    ho.X(:,1,i)=X1n(:,i);  ho.X(:,2,i)=X2n(:,i); 
  end 
end
ho.t=tn; ho.tl=tln; 
p.hopf=ho; % update p.hopf 