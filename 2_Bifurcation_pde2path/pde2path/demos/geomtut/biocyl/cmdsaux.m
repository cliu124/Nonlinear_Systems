%% get branch data a posteriori; here lam1*A: 
dirlist={'l0','l1','l1qr','l2','l2qr','l3qr','l4qr'}; %,'l3qr','l4qr'}; 
nbr=length(dirlist);  bdat=[]; 
for j=1:nbr % loop over branches    
   dir=dirlist{j};  labs=sort(getlabs(dir)); fp=1; incr=1; 
   np=floor(length(labs)/incr)-1; blengths(j)=np; 
   for i=0:np % loop through branch      
    p=loadp(dir,['pt' mat2str(labs(fp+i*incr))]);  %pplot(p,fn); 
    dd=hcylbra2(p,p.u); bdat(j,i+1,:)=dd; 
   end
end
%%
mclf(10); cl={'blue','red','red','m','m','[1 0.6 0]','[0 0.4 0]'}; ax=[-0.8,0.4,-11,10]; 
c1=2; xlab='\lambda_1'; c2=5; ylab='A';  c2=8; ylab='\lambda_1 A'; 
for j=1:nbr; 
  plot(bdat(j,1:blengths(j),c1), bdat(j,1:blengths(j),c2),'color',cl{j},'linewidth',2); 
  hold on; end
xticks([-0.4 0 0.4]); xlabel(xlab); ylabel(ylab); grid on; axis(ax); set(gca,'fontsize',16); 
%% get branch data a posteriori; here lam1*A: 
dirlist={'0l0r','0l1q','0l1qr','0l2qr','0l3qr'}; %,'l3qr','l4qr'}; 
nbr=length(dirlist);  bdat=[]; fn=1; 
for j=1:nbr % loop over branches    
   dir=dirlist{j};  labs=sort(getlabs(dir)); fp=1; incr=1; 
   np=floor(length(labs)/incr)-1; blengths(j)=np; 
   for i=0:np % loop through branch      
    p=loadp(dir,['pt' mat2str(labs(fp+i*incr))]); % pplot(p,fn); 
    dd=hcylbra2(p,p.u); bdat(j,i+1,:)=dd; 
   end
end
%%
mclf(10); cl={'blue','red','red','m','[1 0.6 0]'}; ax=[-0.3,0,-12,0]; 
c1=2; xlab='\lambda_1'; c2=5; ylab='A';  c2=8; ylab='\lambda_1 A'; 
for j=1:nbr; 
  plot(bdat(j,1:blengths(j),c1), bdat(j,1:blengths(j),c2), 'color',cl{j},'linewidth',2); 
  hold on; end
xlabel(xlab); ylabel(ylab); grid on; axis(ax); set(gca,'fontsize',16); 

