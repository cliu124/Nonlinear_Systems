function p=two2lob(p) 
% two2lob: add lobatto points to lin.elements 
gr=p.pdeo.grid; hofem=p.hofem; x=gr.p; 
ne=size(x,2)-1; npoly=p.hofem.femorder*ones(1,ne); 
[xg,c,np]=discr_lob(x,npoly);
hofem.ne=ne; hofem.xe=x; hofem.npoly=npoly; hofem.tri=c; 
gr.p=xg; gr.e(1,2)=np; gr.t=[1:np-1; 2:np]; p.pdeo.grid=gr;  
p.hofem=hofem; p.np=np; p.nu=np*p.nc.neq; 