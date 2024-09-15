function fouplot(p,varargin)
% fouplot: plot of FFT of p.u; currently in 1D and 2D; k-range as aux-arg 

% fouplot(p,{wnr,cmp,kcut}) 
%
% plot FFT of p.u; cmp kcut (dim=1,2) determines a cutoff for plotting in
% Fourier space
[po,tri]=getpte(p); ndim=size(po,1); 
wnr=1; if nargin>1; wnr=varargin{1}; end 
cmp=1; if nargin>2; cmp=varargin{2}; end 
zm=1; if nargin>4; if ~ischar(varargin{4}); zm=varargin{4}; end; end
switch ndim
  case 1; x=po(1,:)'; xmin=min(x); xmax=max(x); nx=2*p.np; 
    xg=linspace(xmin,xmax,nx); 
    kg=-nx/2:1:nx/2-1; % F-vectors (normalized to 2pi) 
    kcut=nx/2; if nargin>3; kcut=varargin{3}; end 
    u1=p.u(1+(cmp-1)*p.np:cmp*p.np); 
    ug=interp1(x,u1,xg); 
    if zm==1; % subtract mean (for better plotting)
     um=sum(sum(ug))/nx; ug=ug-um; fprintf('<u>=%g\n',um); 
    end 
    uf=fftshift(fft(ug))/nx; % normalized FT 
    kf=2*pi/(xmax-xmin); 
    k1=nx/2-kcut+2; k2=nx/2+kcut-1; kvec=kg(k1:k2)*kf;
    figure(wnr); clf; plot(kvec,abs(uf(k1:k2)),'*-'); axis tight; 
    xlabel('k','fontsize',p.plot.fs);ylabel('abs|Fu|','fontsize',p.plot.fs); 
    set(gca,'fontsize',p.plot.fs); 
  case 2; x=po(1,:)'; y=po(2,:)';
    xmin=min(x); xmax=max(x); ymin=min(y); ymax=max(y); nx=4*round(sqrt(p.np)); ny=nx; 
    xg=linspace(xmin,xmax,nx); yg=linspace(ymin,ymax,ny); % "grid" 
    kg=-nx/2:1:nx/2-1; lg=-ny/2:1:ny/2-1; % F-vectors (normalized to 2pi) 
    u1=p.u(1+(cmp-1)*p.np:cmp*p.np); 
    ug=tri2grid(po,tri,u1,xg,yg); % interpol to grid 
    if zm==1; % subtract mean (for better plotting)
      um=sum(sum(ug))/nx/ny; ug=ug-um; fprintf('<u>=%g\n',um); 
    end 
    uf=fftshift(fft2(ug))/nx/ny; % normalized FT 
    kf=2*pi/(xmax-xmin); lf=2*pi/(ymax-ymin); 
    kcut=nx/2*[1,1]; if nargin>3; kcut=varargin{3}; end 
    kvec=(kg(nx/2-kcut(1)+2:nx/2+kcut(1)+1)-0.5)*kf; 
    lvec=(lg(ny/2-kcut(2)+2:ny/2+kcut(2)+1)-0.5)*lf; 
    figure(wnr); clf; 
    pcolor(kvec,lvec,(abs(uf(ny/2-kcut(2)+2:ny/2+kcut(2)+1,nx/2-kcut(1)+2:nx/2+kcut(1)+1)))); 
    axis image; colorbar; grid off; colormap(hot); shading flat;
    set(gca,'fontsize',p.plot.fs)
end
colormap cool;
tit=varargin{end}; if ischar(tit); title(tit,'fontsize',p.plot.fs); end 