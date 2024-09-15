clear all;
close all;

%Estimate region of attraction for
%dx_1=-x_2
%dx_2=x_1+(x_1^2-1)x_2


%We try to find Lyapunov function V as different order of polynomials,
%and let \dot{V}\leq 0 within a
%local region |x|^2\leq \delta^2 using s-procedure, such that 
%As normal, we require V is positive definite, and without loss of
%generality, we can require P>=I. 
sdp_option=sdpsettings('solver','sedumi'); %This solver can be modified as mosek (https://www.mosek.com/), and using edu email, you can get a free license. 

%-----------------------formulate linear matrix inequalities
%define variables to be optimized

degree=16; %polynomial degree

x=sdpvar(2,1);
m = monolist(x,degree/2); %monomials

P = sdpvar(length(m)-1);
V_poly=m(2:end)'*P*m(2:end);

R = sdpvar(length(m));
R_poly=m'*R*m;

I = eye(length(m)-1);

delta2=sdpvar(1,1);

f=([-x(2);x(1)+(x(1)^2-1)*x(2)]);%RHS of nonlinear term. 
dV=jacobian(V_poly, x)*f;
m_p = monolist([x],(degree+2)/2);
Q_V=sdpvar(length(m_p));
Q_dV = sdpvar(length(m_p));
gamma=sdpvar(1);
constraint=[P>=I,...
    coefficients(dV + (delta2-x'*x) * R_poly -m_p'*Q_dV*m_p,m_p)==0, ...
    R>=0,...
    Q_dV<=0];
sol=bisection(constraint, -delta2, sdp_option); %negative sign means maximize delta4. 

% degree_s=4;
% m_s = monolist([x],degree_s/2);
% S = sdpvar(length(m_s));
% S_poly=m_s'*S*m_s;
% 
% beta=sdpvar(1);
% m_con=monolist(x,(degree+degree_s)/2);
% Q_con = sdpvar(length(m_con));
% constraint=[S>=0,...
%     coefficients((value(delta2)-x'*x)-S_poly*(beta-m(2:end)'*value(P)*m(2:end))-m_con'*Q_con*m_con,m_con)==0,...
%     Q_con>=0];
% sol=bisection(constraint,-beta,sdp_option);

x_sym=sym('x',[2,1]);
m_sym = monolist(x_sym,degree/2); %monomials
V_poly=simplify(transpose(m_sym(2:end))*value(P)*m_sym(2:end));
V_fun=matlabFunction(V_poly);

theta_list=linspace(0,2*pi,200);
for theta_ind=1:length(theta_list)
    theta=theta_list(theta_ind);
    x_delta=sqrt(value(delta2))*[cos(theta);sin(theta)];
    V_val_delta(theta_ind)=V_fun(x_delta(1),x_delta(2));
end
[beta,theta_ind]=min(V_val_delta);
theta=theta_list(theta_ind);
x_delta=sqrt(value(delta2))*[cos(theta);sin(theta)];
plot(x_delta(1),x_delta(2),'^k','LineWidth',3); hold on;


x=linspace(-3,3,1000);
y=linspace(-3,3,1000);
[X,Y]=meshgrid(x,y);
V_val=V_fun(X,Y);
V_val(find(V_val>beta))=NaN;
V_val(find(X.^2+Y.^2>value(delta2)))=NaN;
%---------------------
%plot the trajectories from inverse time van der pol oscillator
pcolor(x,y,V_val); shading interp; %colormap(jet);
hold on;

theta_list=linspace(0,2*pi,30);
for theta_ind=1:length(theta_list)
    theta=theta_list(theta_ind);
    x0=0.001*[cos(theta);sin(theta)];
    [t,z]=ode45(@(t,z) [z(2); -z(1)-(z(1)^2-1)*z(2)],[0,20],x0);
    plot(z(:,1),z(:,2),'k','LineWidth',1); hold on;
end

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

xlabel('x_1'); ylabel('x_2');