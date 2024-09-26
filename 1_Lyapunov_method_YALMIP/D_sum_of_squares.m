clear all;
close all;

%Estimate region of attraction for
%dx_1=-x_2
%dx_2=x_1+(x_1^2-1)x_2
%using sum-of-squares

%We try to find Lyapunov function V as different order of polynomials,
%and let \dot{V}\leq 0 within a
%local region |x|^2\leq \delta^2 using s-procedure. 
%As normal, we require V is positive definite, and without loss of
%generality, we can require P>=I. 
sdp_option=sdpsettings('solver','sedumi'); %This solver can be modified as mosek (https://www.mosek.com/), and using edu email, you can get a free license. 

%-----------------------formulate linear matrix inequalities
%define variables to be optimized

degree=4; %polynomial degree of V function

x=sdpvar(2,1);
m = monolist(x,degree/2); %define monomials

P = sdpvar(length(m)-1); %V=m'*P*m, and we only consider the second term to get rid of constant term
V_poly=m(2:end)'*P*m(2:end);

R = sdpvar(length(m)); %This is a SOS multiplier. 
R_poly=m'*R*m;

I = eye(length(m)-1); %identity matrix. 

delta2=sdpvar(1,1); %\dot{V} is negative within |x|^2\leq \delta^2. Her delta2=delta^2.

f=([-x(2);x(1)+(x(1)^2-1)*x(2)]);%RHS of nonlinear term. 

dV=jacobian(V_poly, x)*f;%compute dV/dt. 

m_p = monolist([x],(degree+2)/2);
Q_V=sdpvar(length(m_p));
Q_dV = sdpvar(length(m_p));

constraint=[P>=I,...%V is positive definite
    coefficients(dV + (delta2-x'*x) * R_poly -m_p'*Q_dV*m_p,m_p)==0, ... 
    Q_dV<=0,... %dV is negative definite within a local region |x|^2\leq delta2. 
    R>=0 %make sure R_poly is non-negative multiplier. 
    ]

sol=bisection(constraint, -delta2, sdp_option); %negative sign means maximize delta4. 

%plot the region |x|^2\leq delta^2 as a red line
delta=sqrt(value(delta2));
x_list=linspace(-delta,delta,200);
for x_ind=1:length(x_list)
    x=x_list(x_ind);
    y1(x_ind)=sqrt(delta^2-x^2);
    y2(x_ind)=-sqrt(delta^2-x^2);
end
plot(x_list,y1,'r','Linewidth',3); hold on;
plot(x_list,y2,'r','Linewidth',3); hold on;
pbaspect([1 1 1]);


%get V function as a symbolic function
x_sym=sym('x',[2,1]);
m_sym = monolist(x_sym,degree/2); %monomials
V_poly=simplify(transpose(m_sym(2:end))*value(P)*m_sym(2:end));

V_fun=matlabFunction(V_poly);%convert V_poly to a MATLAB function. 

%get beta=min V within |x|=\delta. 
theta_list=linspace(0,2*pi,200);
for theta_ind=1:length(theta_list)
    theta=theta_list(theta_ind);
    x_delta=delta*[cos(theta);sin(theta)];
    V_val_delta(theta_ind)=V_fun(x_delta(1),x_delta(2));
end
[beta,theta_ind]=min(V_val_delta);

%evaluate V contour for x\in [-3,3] and y\in [-3,3];
x=linspace(-3,3,1000);
y=linspace(-3,3,1000);
[X,Y]=meshgrid(x,y);
V_val=V_fun(X,Y);
V_val(find(V_val>beta))=NaN;% Only consider V\leq \beta
V_val(find(X.^2+Y.^2>delta^2))=NaN; %only consider V that is within the region |x|^2\leq delta^2. 
pcolor(x,y,V_val); shading interp;
hold on;

%plot the trajectories from van der pol oscillator, which reverse the time
%of the given dynamics. 
theta_list=linspace(0,2*pi,30);
for theta_ind=1:length(theta_list)
    theta=theta_list(theta_ind);
    x0=0.001*[cos(theta);sin(theta)]; %initial condition close to origin. 
    [t,z]=ode45(@(t,z) [z(2); -z(1)-(z(1)^2-1)*z(2)],[0,20],x0);  %note that we reverse the time. 
    plot(z(:,1),z(:,2),'k','LineWidth',1); hold on;
end

xlabel('x_1'); ylabel('x_2');