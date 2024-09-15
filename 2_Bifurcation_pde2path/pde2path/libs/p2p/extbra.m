function b1=extbra(b1,b2)
% EXTBRA: append b2 to b1  
try b1=[b1 b2];  % size(b2,1)=size(b1,1)
catch; n1=size(b1,1); n2=size(b2,1); 
  if n1>n2; b2=[b2;zeros(n1-n2,1)]; 
  else; lb1=size(b1,2); b1=[b1; zeros(n2-n1,lb1)]; 
  end
  b1=[b1 b2]; 
end