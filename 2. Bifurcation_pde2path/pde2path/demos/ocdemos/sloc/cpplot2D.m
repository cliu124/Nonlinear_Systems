function cpplot2D(p1,v,su,ps) 
sol=p1.cp;  dia=sldiagn(p1,15); p=p1.oc.s1; np=p.np; nu=p.nu; 
cut=0; sl1=length(sol.t); ts=round(sl1/(su+1)); par=p.u(nu+1:end); 
t=sol.par(1)*sol.t; figure(8); clf; 
iv=[1:ts:su*ts sl1]; ivl=length(iv);  po=getpte(p); 
x1=min(po(1,:)); y1=min(po(2,:)); x2=max(po(1,:)); y2=max(po(2,:)); 
z1a=1e6; z1b=-1e6; z2a=1e6; z2b=-1e6; 
for i=1:ivl 
    u=[sol.u(:,iv(i));par]; p.u=u; %u=real(u); 
    z=slcon(p,u); p.u(np+1:2*np)=z; z1a=min(z1a,min(z)); z1b=max(z1b,max(z)); 
    z2a=min(z2a,min(u(1:np))); z2b=max(z2b,max(u(1:np))); 
    ax(i)=plotsol(p,8,2,ps,'sub',[2,su+1,i]); nola; title(['t=' mat2str(t(iv(i)),3)]); view(v);   
    ax(su+1+i)=plotsol(p,8,1,ps,'sub',[2,su+1,su+1+i]); nola; title([]); view(v); 
end
if ps==2; 
 for i=1:su+1; colormap(ax(i),parula); caxis(ax(i),[z1a, z1b]); end 
 for i=su+2:2*su+1; caxis(ax(i),[z2a, z2b]); end 
else
  ax1=[x1 x2 y1 y2 z1a z1b];   ax2=[x1 x2 y1 y2 z2a z2b];   
  for i=1:su+1; colormap(ax(i),hot); axis(ax(i),ax1); caxis(ax(i),[z1a, z1b])
  end 
  for i=su+2:2*su+2; axis(ax(i),ax2); colormap(ax(i),parula); caxis(ax(i),[z2a, z2b])
end 
end