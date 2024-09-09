%% illustrate full spectrum!
p=loadp('FSS','pt0');
p.nc.neig=p.nu; r=resi(p,p.u); Gu=getGupde(p,p.u,r); 
[ineg,muv]=spcalc(Gu,p);
%% get the defects of CSS
p=loadp('FSS','pt2'); [Psi,muv,d,t1]=getPsi(p);
p=loadp('FSS','pt6'); [Psi,muv,d,t1]=getPsi(p);
%% floqps works reasonably 
[muv1, muv2,ind,h]=floqpsap('h2','pt4'); fprintf('d(u_H)=%i\n',ind); 