% convenience function to load CP previously saved via savecp
function [fn,alv,vv,tlv,sol,tv,uv,opt]=loadcp(fname)
global s0 s1 u1 Psi par xi um1 um2 sig;
p=load(fname); fn=p.fn; 
s1=loadp(fn.sd1,fn.sp1); s0=loadp(fn.sd0,fn.sp0); par=getaux(s1); 
n=s1.nu; u1=s1.u(1:n); opt=p.opt; 
Psi=p.Psi; xi=p.xi; alv=p.alv; vv=p.vv; tlv=p.tlv; 
sol=p.sol; 
tv=p.tv; uv=p.uv; %um1=p.um1; um2=p.um2; 