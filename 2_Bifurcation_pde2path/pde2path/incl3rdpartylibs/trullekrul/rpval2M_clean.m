function out = rpval2M_clean(vin)
% takes out doubles
vin = sort(vin,2);
I = and([true(size(vin,1),1) vin(:,1:end-1)~=vin(:,2:end)],vin~=0);
vin(not(I)) = 0;
vin = sort(vin,2);
nN = sum(vin~=0,2);
out = vin(:,size(vin,2):-1:size(vin,2)-max(nN)+1);