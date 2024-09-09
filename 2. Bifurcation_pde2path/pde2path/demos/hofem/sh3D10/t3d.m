[ne,ng,po,c,efl,gfl] = cube10 (0);
%%
for i=1:ne
    x=po(c(i,:),1); y=po(c(i,:),2); z=po(c(i,:),3); f=0*x; 
    volume = t10vis (x,y,z,f,1); view(3); axis([-1 3 -1 1 -1 1]); hold on
    pause 
end 
%%
[ne,ng,po,c,efl,gfl] = cube10b(2,0); ne
%%
figure(1); clf; 
po=p.pdeo.grid.p'; c=p.tri; 
for i=1:p.nt
    x=po(c(i,:),1); y=po(c(i,:),2); z=po(c(i,:),3); f=0*x; i
    volume = t10vis (x,y,z,f,1); view(3); %axis([-1 3 -1 1 -1 1]); 
    hold on
    pause; %(0.01);  
end 