%Baleanu, D., Jajarmi, A., & Hajipour, M. (2018).
%On the nonlinear dynamical systems within the generalized fractional derivatives with Mittag-Leffler kernel
%Nonlinear dynamics, 94(1), 397-414.
% The numerical coded in this m file is given in eqs.(33)-(34) in the above
% research paper
clc; clear; close all;
% Initial Conditions
t(1)=0; x(1)=3.5;y(1)=1.5;
% Inputs
alpha=0.892;ABC=1;h= 0.01; tfinal=100; t=t(1):h:tfinal;
% Parameters of the system
g=0.011;theta=0.0251;eta=3.36;C=0.51;varep=0.001;
% Number of Iterations
N=ceil((tfinal-t(1))/h);
% The given fractional-order dynamical system under the ABC operator
f1 = @(t,x,y) -g*x-theta*y;
f2= @(t,x,y) eta*x+C*y-varep*y.^3;
% ABC Fractional Euler Method starts
tic;
for n=1:N
j=1:n;
x(n+1)=x(1)+((1-alpha)*f1(t(n),x(n),y(n)))/(ABC)+...
((alpha.*h^alpha)./(ABC.*gamma(alpha+1))).*...
sum(((n-j+1).^alpha-(n-j).^alpha).*f1(t(j),x(j),y(j)));

y(n+1)=y(1)+((1-alpha)*f2(t(n),x(n),y(n)))/(ABC)+...
((alpha.*h^alpha)./(ABC.*gamma(alpha+1))).*...
sum(((n-j+1).^alpha-(n-j).^alpha).*f2(t(j),x(j),y(j)));

t(n+1)=t(n)+h;
end
toc;
figure(1)
plot(t,x,'k-')
hold on
plot(t,y,'r-'),xlabel('time'),ylabel('Approximation for x(t) and y(t)'),legend('x_t','y_t')
figure(2)
plot(x,y),xlabel('x-values'),ylabel('y-values'),legend('phase plane')