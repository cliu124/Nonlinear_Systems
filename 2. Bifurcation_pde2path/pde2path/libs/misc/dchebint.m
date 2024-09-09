function I=dchebint(p,u)
% dchebint: integrate u (in ChebFourier-discr.) over disk
% p.na=#points in angle, p.r=radial Chebychev points (in (0,1])
% disk-radius assumed in p.lx
uu=reshape(u,p.na,p.nr); I=0; 
for i=1:p.na;  % loop over angles 
   I=I-trapz(p.lx*p.r,p.lx*p.r'.*uu(i,:)); % wrong order of r, hence -
end 
I=2*pi*I/p.na; 