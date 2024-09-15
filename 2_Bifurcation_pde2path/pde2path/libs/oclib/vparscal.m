function s=vparscal(v1,par1,v2,par2,c1,c2)
% vparscal: scaling of parameters for norm, used in isc  
s=c1*(sum(v1(:).*v2(:))+c2*par1(1)*par2(1))+(1-c1)*par1(2)*par2(2);