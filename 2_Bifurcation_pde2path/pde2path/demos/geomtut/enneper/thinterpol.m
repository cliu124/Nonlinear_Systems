function p=thinterpol(p,idold,tho) % interpolate theta to new bdry nodes 
idnew=p.idx; % current bdry ids 
nn=setdiff(idnew,idold,'stable'); % new nodes 
nl=length(nn); %size(idnew), size(idold), nl, pause 
if nl>0; size(idold), size(idnew), end
thn=tho;  X=p.X; 
for i=1:nl % loop over new bd points and augment thn by new th-vals 
      na=nn(i); % new bd point index         
      Xd=dot(X(idold,:)-X(na,:),X(idold,:)-X(na,:),2); 
      [Xds,sid]=sort(Xd,'ascend');          
      th1=tho(idold(sid(1))); th2=tho(idold(sid(2))); 
      newt=0.5*(th1+th2); 
      if th1*th2<0; 
         if th1>pi/2; 
             if newt>0; newt=pi-newt; else; newt=newt+pi; end
         else 
         end
      end    
      thn=[thn; newt]; 
end 
p.th=thn;  