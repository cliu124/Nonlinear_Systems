function  [K,M,Kx]=assem1Dq(p) 
% assem1Dq: assemble K,M,Kx for 1D quadratic elements
% after Pozrikidis. elem-end-points xe in p.hofem.xe
% hofem and grid not modified; preparation by lin2quad 
xe=p.hofem.xe; ne=size(xe,2)-1; 
for l=1:ne;   h(l)=xe(l+1)-xe(l); end % elem size 
ng=2*ne+1; 
% diagonals of pentadiagonal matrices
ap=zeros(ng,1); bp=ap; cp=ap; dp=ap; ep=ap; 
apm=ap; bpm=ap; cpm=ap; dpm=ap; epm=ap; 
apc=ap; bpc=ap; cpc=ap; dpc=ap; epc=ap; 
% loop over all elements
for l=1:ne    
  cf=1.0/(3.0*h(l));
  A11=7.0*cf;  A12=-8.0*cf; A13=cf;
  A21=A12;     A22=16.0*cf; A23=A12;
  A31=A13;     A32=A23;     A33=A11;

  cf=h(l)/30.0;
  B11=4.0*cf; B12= 2.0*cf;  B13=-cf;
  B21=B12;    B22=16.0*cf;  B23=B12;
  B31=B13;    B32=B23;      B33=B11;
   
  C11=-1/2; C12= 2/3;  C13=-1/6; 
  C21=-C12;  C22=0;  C23=2/3;
  C31=-C13;  C32=-C23;  C33=-C11;

  cl1=2*l-1; cl2=2*l; cl3=2*l+1;
  ap(cl1)=ap(cl1)+A11; bp(cl1)=bp(cl1)+A12; cp(cl1)=cp(cl1)+A13;  
  dp(cl2)=dp(cl2)+A21; ap(cl2)=ap(cl2)+A22; bp(cl2)=bp(cl2)+A23;
  ep(cl3)=ep(cl3)+A31; dp(cl3)=dp(cl3)+A32; ap(cl3)=ap(cl3)+A33;
  
  apm(cl1)=apm(cl1)+B11; bpm(cl1)=bpm(cl1)+B12; cpm(cl1)=cpm(cl1)+B13;  
  dpm(cl2)=dpm(cl2)+B21; apm(cl2)=apm(cl2)+B22; bpm(cl2)=bpm(cl2)+B23;
  epm(cl3)=epm(cl3)+B31; dpm(cl3)=dpm(cl3)+B32; apm(cl3)=apm(cl3)+B33;
  
  apc(cl1)=apc(cl1)+C11; bpc(cl1)=bpc(cl1)+C12; cpc(cl1)=cpc(cl1)+C13;  
  dpc(cl2)=dpc(cl2)+C21; apc(cl2)=apc(cl2)+C22; bpc(cl2)=bpc(cl2)+C23;
  epc(cl3)=epc(cl3)+C31; dpc(cl3)=dpc(cl3)+C32; apc(cl3)=apc(cl3)+C33;  
end
K=spdiags([cp,bp,ap,dp,ep],-2:2,ng,ng); K=K';  
k=1; M=spdiags([cpm,bpm,apm,dpm,epm],-2:2,ng,ng); M=M'/k;  
Kx=spdiags([cpc,bpc,apc,dpc,epc],-2:2,ng,ng); Kx=Kx';  
