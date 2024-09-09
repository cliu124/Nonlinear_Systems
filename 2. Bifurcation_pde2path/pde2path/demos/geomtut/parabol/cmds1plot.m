%% plotting of diagnostics for parabol problem 
%% raw and loglog, e_0 only useful for 5,6 
mclf(4); c=5;  xti=[100 500 2500]; yti='auto'; %ev2=ev3; 
switch c    
    case 4; ylab='||G||_\infty';  
    case 5; ylab='||z-z_a||_\infty';   
    case 6; ylab='||z-z_a||_2';     
    case 11; ylab='\delta';     
end 
si=2; e0=ev0(:,si:end); e1=ev1(:,si:end); e2=ev2(:,si:end); 
plot(e0(1,:),e0(c,:),'-*'); hold on; plot(e1(1,:),e1(c,:),'-+'); 
plot(e2(1,:),e2(c,:),'-o'); xticks(xti); % raw data 
xlabel('n');  ylabel(ylab); legend('S_0','S_1','S_2'); 
set(gca,'fontsize',14); %axis tight; grid on; 
%
mclf(5); loglog(e0(1,:),e0(c,:),'-*'); hold on; % loglog-plots 
loglog(e1(1,:),e1(c,:),'-x'); loglog(e2(1,:),e2(c,:),'-o'); 
logn0=log(e0(1,:)); loge0=log(e0(c,:)); yti=0.00005*[1 10 100]; 
pf=polyfit(logn0, loge0, 1); al0=pf(1);  % fit conv.rate 
logn1=log(e1(1,:)); loge1=log(e1(c,:)); 
pf=polyfit(logn1, loge1, 1); al1=pf(1);  
logn2=log(e2(1,:)); loge2=log(e2(c,:)); 
pf2=polyfit(logn2, loge2, 1); al2=pf2(1); 
%y=polyval(pf2,logn2); k2=exp(pf2(2)); loglog(ev2(1,:),exp(y)); % plot fit
legend(['S_0, \alpha=' mat2str(al0,3)],['S_1, \alpha=' mat2str(al1,3)],...
    ['S_2, \alpha=' mat2str(al2,3)])
xlabel('n'); ylabel(ylab); 
set(gca,'fontsize',14); axis tight; grid on; xticks(xti); yticks(yti); 
%% some norms 
mclf(4); c=5;  xti=[50 500 2500]; yti='auto'; ev1=ev2; 
switch c    
    case 4; ylab='||G||_\infty';  
    case 5; ylab='||z-z_a||_\infty';   
    case 6; ylab='||z-z_a||_2';     
    case 11; ylab='\delta';     
end 
mclf(5); loglog(e0(1,:),e0(c,:),'-*'); hold on; loglog(e1(1,:),e1(c,:),'-x'); 
logn0=log(e0(1,:)); loge0=log(e0(c,:)); 
yti=round([exp(loge0(end-1)),exp(loge0(3))],1,'significant'); 
pf=polyfit(logn0, loge0, 1); al0=pf(1);  
%y=polyval(pf,logn0); k0=exp(pf(2)); loglog(ev0(1,:),exp(y)); 
logn1=log(e1(1,:)); loge1=log(e1(c,:)); 
pf=polyfit(logn1, loge1, 1); al1=pf(1);  
%y=polyval(pf,logn0); k0=exp(pf(2)); loglog(ev0(1,:),exp(y)); 
legend(['e_0, \alpha=' mat2str(al0,3)],['e_1, \alpha=' mat2str(al1,3)]); 
xlabel('n'); ylabel(ylab); set(gca,'fontsize',14); axis tight; grid on; 
xticks(xti); yticks(yti); 
%% quotient e2/e0
q=ev2(c,:)./ev1(c,:); mclf(6); semilogx(ev2(1,:),q); axis tight; grid on; 
%% other stuff, only e1,e2 useful 
mclf(4); c=11;  xti=[50 500 2500]; yti='auto'; 
switch c
    case 2; ylab='res';    case 4; ylab='|| G||_\infty'; 
    case 5; ylab='|| z-z_a||_\infty';  case 6; ylab='|| z-z_a||_2';  
    case 7; ylab='|| H-H_a||_\infty';  case 8; ylab='|| H-H_a||_2';  
    case 9; ylab='|| H-H_a||_{\infty,weak}'; case 10; ylab='|| K_f-1||_2';         
    case 11; ylab='\delta'; case 15; ylab='h';  case 16; ylab='<node valence>';   
end 
plot(ev1(1,:),ev1(c,:),'-x'); hold on; plot(ev2(1,:),ev2(c,:),'-o'); 
xlabel('n');  ylabel(ylab); legend('e_1','e_2'); 
set(gca,'fontsize',14); axis tight; grid on; xticks(xti); yticks(yti); 
%%
mclf(5); loglog(ev1(1,:),ev1(c,:),'-x'); hold on; loglog(ev2(1,:),ev2(c,:),'-o'); 
logn1=log(ev1(1,:)); loge1=log(ev1(c,:)); 
yfac=1e5; yti=round(yfac*[exp(loge1(end-1)),exp(loge1(2))])/yfac; 
pf=polyfit(logn1, loge1, 1); al1=pf(1);  
logn2=log(ev2(1,:)); loge2=log(ev2(c,:)); 
pf2=polyfit(logn2, loge2, 1); al2=pf2(1); 
legend(['e_1, \alpha=' mat2str(al1,3)],['e_2, \alpha=' mat2str(al2,3)]); 
xlabel('n'); ylabel(ylab); 
set(gca,'fontsize',14); axis tight; grid on; xticks(xti); yticks(yti); 