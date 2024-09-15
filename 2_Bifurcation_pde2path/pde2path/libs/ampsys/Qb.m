function qqd=Qb(idx,qb,qbc,qbcc,sk)
% Qb: helper function for ampsys
i1=idx(1);i2=idx(2);
if i1<=sk && i2<=sk
    qqd=[i1,i2,qb];
end
if i1<=sk && i2>sk
   qqd=[i1,i2,qbc]; 
end
if i1>sk && i2>sk
   qqd=[i1,i2,qbcc]; 
end
