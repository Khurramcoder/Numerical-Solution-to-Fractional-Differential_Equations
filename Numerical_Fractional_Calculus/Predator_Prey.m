% The Fractional Forward Euler's Methodfor a System of ODEs
% Predator Prey System
%Dr. Sania Qureshi

clc; clear; close all;

%Inputs

h=0.01; t(1)=0;x(1)=1; y(1)=3; tfinal=10; t=t(1):h:tfinal; N=ceil((tfinal-t(1))/h);
% Initial Conditions

% Constant fractional order
alpha=0.95; beta=0.90;


%Given below is the Fractional Order System of ODEs(Predator Prey System)
f1= @(t,x,y)x.*(4-y);
f2= @(t,x,y)y.*(x-5);

% The Fractional Forward Euler Method for Fractional Order System of ODEs(Predator Prey System)

for n=1:N
    j=1:n;
    t(n+1)=t(n)+h;
    x(n+1)=x(1)+((h^alpha)/(gamma(alpha+1))).*sum(((n-j+1).^(alpha)-(n-j).^(alpha)).*f1(t(j),x(j),y(j)));
    y(n+1)=y(1)+((h^beta)/(gamma(beta+1))).*sum(((n-j+1).^(beta)-(n-j).^(beta)).*f2(t(j),x(j),y(j)));
end
figure(1)
plot(t,x)
figure(2)
plot(t,y)
figure(3)
plot(x,y)