bw=3; beltol=1e-6; belimax=5; % border-width, bel-parameters 
AMG=0; % set AMG=1 if ilupack is available 
if ~AMG; p=setbel(p,bw,beltol,belimax,@lss); % use BEL without ilupack 
else p=setbel(p,bw,beltol,belimax,@lssAMG);  end 