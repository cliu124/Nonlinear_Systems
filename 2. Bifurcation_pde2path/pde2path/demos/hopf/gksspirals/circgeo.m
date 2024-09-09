function geo=circgeo(r,nx)
td=linspace(0,2*pi,nx); geo=polygong(r*cos(td),r*sin(td));  