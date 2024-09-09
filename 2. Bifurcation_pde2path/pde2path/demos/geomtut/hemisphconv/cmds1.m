%% hemisphere conv-test 
keep pphome; hucl; global p2pglob; p2pglob.edc='k'; p2pglob.cut=0; p2pglob.showN=0; 
p2pglob.cb=1; p2pglob.axlab=1; p2pglob.showbd=2; 
par=-1; p=[]; p=sphereinit(par,1,1); pplot(p); 
%% loop over different refinements
ev=[]; p.sw.verb=2; nref=5; 
for i=1:nref; p=sphereinit(par,1,i); ev=[ev eplot(p)];  end
%% plot, see also cmds1plot for more sophisticated plotting 
mclf(4); c=3;  xti=[50 500 5000];  yti='auto'; 
switch c
    case 2; ylab='res'; 
    case 3; ylab='|| H+1||_\infty';  case 4; ylab='|| H_f+1||_\infty'; 
    case 5; ylab='|| K-1||_\infty';  case 6; ylab='|| K_f-1||_\infty';  
    case 7; ylab='|| H+1||_2';     case 8; ylab='|| H_f+1||_2';  
    case 9; ylab='|| K-1||_2';    case 10; ylab='|| K_f-1||_2';      
    case 12; ylab='|| |X|-1||_\infty';  case 13; ylab='|| |X|-1 ||_2';  
    case 15; ylab='|| |X_f|-1||_\infty'; case 16; ylab='|| |X_f|-1 ||_2'; 
    case 17; ylab='\delta'; case 20; ylab='h';  case 22; ylab='<node valence>';   
end 
plot(ev(1,:),ev(c,:),'-*'); xlabel('n');  ylabel(ylab); % raw 
set(gca,'fontsize',14); axis tight; grid on; 
mclf(5); loglog(ev(1,:),ev(c,:),'-*'); hold on % loglog 
xlabel('n'); ylabel(ylab); 
set(gca,'fontsize',14); axis tight; grid on; xticks(xti); 
logn1=log(ev(1,:)); loge1=log(ev(c,:)); 
pf=polyfit(logn1, loge1, 1); y=polyval(pf,logn1); al1=pf(1); k1=exp(pf(2)); 
loglog(ev(1,:),exp(y)); yticks(yti); als=mat2str(al1,3); 
legend(['\alpha=' als]); 