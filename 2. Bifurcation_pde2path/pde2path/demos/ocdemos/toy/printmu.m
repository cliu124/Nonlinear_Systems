function printmu(p) % auxx to print multipliers 
try; mu=p.oc.muv; catch; mu=p; end; n1=length(mu); 
fprintf(['leading multipliers\n' repmat('%g ',1,n1) '\n'], mu(1:end)); 