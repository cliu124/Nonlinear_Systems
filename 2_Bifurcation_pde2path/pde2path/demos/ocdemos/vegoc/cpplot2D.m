function cpplot2D(p1,v,su,ps,pfak,fn) 
load('vegcm.asc'); load('watcm.asc'); sol=p1.cp; p=p1.oc.s1; np=p.np; nu=p.nu; 
vegdiagn(p1,15,pfak,fn); sl1=length(sol.t); ts=round(sl1/(su+1)); par=p.u(nu+1:end); 
t=sol.par(1)*sol.t; figure(8); clf; 
iv=[1:ts:su*ts sl1]; ivl=length(iv);  po=getpte(p); 
x1=min(po(1,:)); y1=min(po(2,:)); x2=max(po(1,:)); y2=max(po(2,:)); 
z1a=1e6; z1b=-1e6; z2a=1e6; z2b=-1e6; z3a=1e6; z3b=-1e6; su=su+1; 
for i=1:ivl     
    u=[real(sol.u(:,iv(i)));par]; p.u=u; z=efu(p); 
    p.u(2*np+1:3*np)=z; z1a=min(z1a,min(z)); z1b=max(z1b,max(z)); 
    z2a=min(z2a,min(u(1:np))); z2b=max(z2b,max(u(1:np))); 
    z3a=min(z3a,min(u(np+1:2*np))); z3b=max(z3b,max(u(np+1:2*np))); 
    ax(i)=plotsola(p,8,3,ps,'sub',[3,su,i]); nolti; title(['t=' mat2str(t(iv(i)),3)]); view(v);   
    ax(su+i)=plotsola(p,8,1,ps,'sub',[3,su,su+i]); nolti; title([]); view(v); 
    ax(2*su+i)=plotsola(p,8,2,ps,'sub',[3,su,2*su+i]); nolti; title([]); view(v); 
end
if ps==2; 
 for i=1:su+1; colormap(ax(i),vegcm); caxis(ax(i),[z1a, z1b]); end 
 for i=su+2:2*su+1; caxis(ax(i),[z2a, z2b]); end 
 for i=2*su+3:3*su+1; caxis(ax(i),[z3a, z3b]); end 
else
  ax1=[x1 x2 y1 y2 z1a z1b];   ax2=[x1 x2 y1 y2 z2a z2b]; ax3=[x1 x2 y1 y2 z3a z3b];   
  for i=1:su; colormap(ax(i),hot); axis(ax(i),ax1); caxis(ax(i),[z1a, z1b]); 
      zticks(ax(i),[ceil(z1a) floor(z1b)]); 
  end 
  for i=su+1:2*su; axis(ax(i),ax2); colormap(ax(i),vegcm); caxis(ax(i),[z2a, z2b]);
      zticks(ax(i),[ceil(z2a) floor(z2b)]); 
  end 
  %watcm(:,1)=watcm(:,1)/4; 
  for i=2*su+1:3*su; axis(ax(i),ax3); colormap(ax(i),watcm); caxis(ax(i),[z3a, z3b]);
      zticks(ax(i),[ceil(z3a) floor(z3b)]); 
  end 
end