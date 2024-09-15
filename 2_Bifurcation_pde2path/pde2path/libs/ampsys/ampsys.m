function [Q,C,c1,phi]=ampsys(p)
% ampsys: compute coefficients for Landau system associated to Turing-Bif
%
% struct p contains:
% the (critical) wave vector lattice in p.k
% 
% for scalar (SH-type) equations: the quadratic and cubic coeff. in p.c2
%    and p.c3; if p.c2 is not given, then aecoeff assumes a function c2=c2(k) 
%    in the current directory. Additionally, l=L(k) must give the Fourier symbol 
%    of the linear part
% for RD system: the diffusion matrix p.D, the spatially homog. sol. p.uh, 
%    the bifurcation param p.bifpar and p.uhp=uh(par(bifpar)+del)  
%    the nonlinearity f=f(u,p) must be in the current dir. 
% optional: p.eqnr,  list of eqns for which to compute coeff., 
%    default  p.eqnr=1
% switches: 
%   p.sb: if 0, then purely numerical; if 1, then symbolic (SH-case)
%   p.cons: if 0, then solve for O(eps^2) terms at regular wave-vec (default)
%           if 1, then drop all quadr.terms except those at k_1,...,k_m 
%   p.qs: function handle for Fourier symbol of quadr. terms, e.g.   KS
zerocut=1e-6; Co=[]; Qo=[]; 
try; eqnr=p.eqnr; catch eqnr=1; end 
try; cons=p.cons; catch cons=0; end 
try; sb=p.sb; catch; sb=0; end  % sb=0 means no symbolic expressions 
qs=0; if isfield(p,'qs'); qs=1; end
if isfield(p,'D'); cmp=size(p.D,1);D=p.D; uh=p.uh; % check SH or RD 
else; cmp=1; c3=p.c3; c1=0; % return dummy for c0
end
sk=size(p.k,2); k=[p.k,-p.k]; % the wave-vectors
fprintf('k=');     
for i1=1:size(k,1); % rows 
    fprintf('\n');    
    for i2=1:sk; fprintf('%s',printcon(k(i1,i2))); end
end
fprintf('\n'); 
if cmp>=2  % if 2 comp. system, build J, L, phi, B, T  
 u=sym('u',[cmp 1]); Jf=jacobian(f(u,p),u);  Js=subs(Jf,u,uh); 
 kc=norm(k(:,1)); LL=Js-D*norm(kc)^2; LL=double(LL); [ev,ew]=eig(LL); 
 w=zeros(cmp,1);  for i=1:cmp; w(i)=ew(i,i); end; [~,idx]=sort(abs(real(w)));
 disp(['zero eigenvalue: ',num2str(ew(idx(1),idx(1)))]); 
 p.par(p.bifpar)=p.par(p.bifpar)+p.del; Jfp=jacobian(f(u,p),u); % Jfp, u, p.uhp, pause 
 Jsp=subs(Jfp,u,p.uhp);
 LLp=Jsp-D*norm(kc)^2; LLp=double(LLp);  p.par(p.bifpar)=p.par(p.bifpar)-p.del;
 [~,ewp]=eig(LLp); wp=w;  for i=1:cmp; wp(i)=ewp(i,i); end; [~,idxp]=sort(abs(real(wp))); 
 c1=(ewp(idxp(1),idxp(1))-ew(idx(1),idx(1)))/p.del; 
 fprintf('c1=%g\n',c1); 
 phi=ev(:,idx(1)); if phi(1)~=0; phi=phi/phi(1); else phi=phi/norm(phi); end; 
 LT=conj(transpose(LL)); [ev,ew]=eig(LT);
 w=zeros(cmp,1); for i=1:cmp; w(i)=ew(i,i); end
 [~,idx]=sort(abs(real(w))); psi=ev(:,idx(1));
 psi=psi/(phi'*psi); b=asB(phi,phi,p,cmp); bc=asB(conj(phi),phi,p,cmp); 
 bcc=asB(conj(phi),conj(phi),p,cmp); 
 t=psi'*asT(phi,phi,phi,p,cmp); tc=psi'*asT(conj(phi),phi,phi,p,cmp);
 tcc=psi'*asT(conj(phi),conj(phi),phi,p,cmp); tccc=psi'*asT(conj(phi),conj(phi),conj(phi),p,cmp);
end
for je=eqnr  % loop over eq-nr: compute coeff. in \pa_t A_je=...
fprintf('eq-nr %i:\n', je);  
m=size(k,2); qd=[]; dq=1; Q=[]; % find quadratic terms
for i1=1:m
    for i2=1:m;        
         if k(:,je)==k(:,i1)+k(:,i2) % cannot remove this term
            dq=0;
            if cmp>=2; Q=[Q;Qb(sort([i1 i2]),psi'*b,psi'*bc,psi'*bcc,sk)];
            else
              if qs; Q=[Q;sort([i1 i2]) p.c2*p.qs(k(:,i1)+k(:,i2))];
              else Q=[Q;sort([i1 i2]) p.c2];   end                
            end
         else  qd=[qd; sort([i1 i2]) 1]; % can remove this term 
         end
    end
end
if dq==0 % cannot remove all quad terms, inform user 
    Q=myfilter(Q);
    if je==1; disp('Cannot remove quadratic terms from the residual, since');
    for i=1:size(Q,1)
    disp(['k',num2str(double(Q(i,1))),'+k',num2str(double(Q(i,2))),'=k1']);
    end
    if cons; qd=[]; end;   % if consistent expansion, then remove qd 
    end
end
if ~cons; dq=1; end   % solve for quadr.terms where possible
if dq==1;  
  qd=myfilter(qd); % can remove all quad terms
  qqd=zeros(size(qd,1),size(qd,2)+cmp-1);
  if sb==1; qqd=sym(qqd); end
  for i=1:size(qd,1)
    if cmp>=2;  qqd(i,:)=Qb(qd(i,1:2),b'*qd(i,3),bc'*qd(i,3),bcc'*qd(i,3),sk);
    else
      if qs; qqd(i,:)= [qd(i,1:2), p.c2*p.qs(k(:,qd(i,1))+k(:,qd(i,2)))*qd(i,3)]; 
      else qqd(i,:)= [qd(i,1:2),p.c2*qd(i,3)]; end            
    end
  end
  qd=qqd;
  for i=1:size(qd,1) % apply L^-1 (if L~=0) to the quad terms and to the ansatz
    i1=qd(i,1);i2=qd(i,2);
    if cmp==1; l=L(k(:,i1)+k(:,i2),p);  
    else; l=Js-D*norm(k(:,i1)+k(:,i2))^2; end
    if abs(det(l))>zerocut; d=-l\qd(i,3:end)'; qd(i,3:end)=d'; end; 
  end
  %qd, pause
  % calc Q(oldansatz,newansatz) and find terms corresponding to k1
  er=[];
  for i=1:size(k,2)
   for j=1:size(qd,1)
    kk=k(:,i)+k(:,qd(j,1))+k(:,qd(j,2));
    if kk==k(:,je)
        s=sort([i qd(j,1) qd(j,2)]);
        if cmp>=2
            if i<=sk; er=[er; s 2*psi'*asB(phi,qd(j,3:end)',p,cmp)];
            else; er=[er; s 2*psi'*asB(conj(phi),qd(j,3:end)',p,cmp)]; end
        else
          if qs;  er=[er; s 2*p.c2*p.qs(k(:,je))*qd(j,3)]; 
          else  er=[er; s 2*p.c2*qd(j,3)];  end            
        end
    end
   end 
 end
qd=myfilter(er);
end  % if dq==1, 
% calc C(oldansatz,oldansatz,oldansatz) and find terms corresponding to k1
m=size(k,2); er=[];
for i1=1:m
  for i2=1:m
     for i3=1:m
        kk=k(:,i1)+k(:,i2)+k(:,i3);
        if kk==k(:,je)
           if cmp>=2;  er=[er;Ct(sort([i1 i2 i3]),t,tc,tcc,tccc,sk)];
           else; er=[er;sort([i1 i2 i3]) c3]; end
        end
      end
  end
end
cb=myfilter(er);
% compare qd and cd and put together
if dq==1
  for i=1:size(qd,1)
    z=qd(i,4);ix=[];
    for j=1:size(cb,1)      
        if qd(i,1:3)==cb(j,1:3); z=z+cb(j,4); ix=[ix,j];  end 
    end
    if isempty(ix)==0; cb(ix,:)=[]; end    
    qd(i,4)=z;
  end
end
if cons==1; C=cb; else; C=[qd;cb]; end % HU 
if cmp==1;phi=1;end
% change something like C=(1 1 2 5) into C=(1 1 -1 5) 
for i=1:size(C,1)
    for j=1:size(C,2)-1
        if C(i,j)>sk, C(i,j)=-(C(i,j)-sk); end
    end
end
% the same for Q
if Q~=0
    for i=1:size(Q,1)
      for j=1:size(Q,2)-1
        if Q(i,j)>sk; Q(i,j)=-(Q(i,j)-sk); end
      end
    end 
end
if size(Q,2)>1; fprintf('Q=\n'); 
   for i1=1:size(Q,1)
    if sb==0; fprintf('%i %i %g\n', Q(i1,1), Q(i1,2), Q(i1,3));
    else fprintf('%i %i %s\n', Q(i1,1), Q(i1,2), Q(i1,3)); end 
   end
end
fprintf('C=\n');     
for i1=1:size(C,1); 
    if sb==0; fprintf('%i %i %i %4.3f\n', C(i1,1), C(i1,2), C(i1,3), C(i1,4)); 
    else fprintf('%i %i %i %s\n', C(i1,1), C(i1,2), C(i1,3), C(i1,4)); end 
end
Co=[Co C]; Qo=[Qo Q]; 
end
C=Co; Q=Qo; 
