function s=Ct(idx,t,tc,tcc,tccc,sk)
% Ct: helper function in ampsys
i1=idx(1);i2=idx(2);i3=idx(3);
if i1<=sk && i2<=sk && i3<=sk
   s=[idx,t];  
end
if (i1<=sk && i2<=sk && i3>sk) || (i1<=sk && i2>sk && i3<=sk) || (i1>sk && i2<=sk && i3<=sk)  
   s=[idx,tc];  
end
if (i1<=sk && i2>sk && i3>sk) || (i1>sk && i2>sk && i3<=sk) || (i1>sk && i2<=sk && i3>sk)  
   s=[idx,tcc];  
end
if (i1>sk && i2>sk && i3>sk)
   s=[idx,tccc];  
end