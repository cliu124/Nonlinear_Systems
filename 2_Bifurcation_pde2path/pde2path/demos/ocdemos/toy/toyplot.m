function toyplot(poc)
for k=1:4; figure(k);clf; plot(poc.cp.par(1)*poc.cp.t,poc.cp.u(k,:)); 
    xlabel('t'); ylabel(['u_' mat2str(k)]); axis tight; set(gca,'fontsize',14); end