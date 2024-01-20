% Qureshi,S., Atangana,A., & Shaikh,A. A.(2019)
% Strange Chaotic attractors under fractal-fractional operators using newly
% proposed numerical methods.
% The European Physical Journal Plus, 134(10),523.
% This numerical scheme is for Fractal-Fractional Caputo Fractional
% Differential Operator
%This scheme is proposed by Prof. Dr. Atangana.


clc;clear;close all;
h=0.01;
t(1)=0.1;
x(1)=-8.1;y(1)=-2;z(1)=-7.1;
alpha=0.99;
tau=0.95;
AB=1-alpha+alpha/gamma(alpha);
tfinal=200;
t=t(1):h:tfinal;
N=ceil((tfinal-t(1))/h);
lambda=-5;epsilon=-6;c=5;d=-1;a=7.5;b=1;
%Nonlinear System of ODEs
f1=@(t,x,y,z)a*(y-x)+b*y.*(z.^2);
f2=@(t,x,y,z)c*x+d*x.*(z.^2);
f3=@(t,x,y,z)lambda*z+epsilon*abs(x);
% Algorithm Starts
for n=1:N
    j=2:n;
    x(n+1)=x(1)+(tau*(h^alpha)./gamma(alpha+2)).*sum(((n+1-j).^alpha.*...
        (n-j+2+alpha)-(n-j).^alpha.*(n-j+2+2*alpha)).*...
        f1(t(j),x(j),y(j),z(j)).*t(j).^(tau-1)-((n+1-j).^(alpha+1)-(n-j).^alpha.*...
        (n-j+1+alpha)).*f1(t(j-1),x(j-1),y(j-1),z(j-1)).*t(j-1).^(tau-1));
    
    y(n+1)=y(1)+(tau*(h^alpha)./gamma(alpha+2)).*sum(((n+1-j).^alpha.*...
        (n-j+2+alpha)-(n-j).^alpha.*(n-j+2+2*alpha)).*...
        f2(t(j),x(j),y(j),z(j)).*t(j).^(tau-1)-((n+1-j).^(alpha+1)-(n-j).^alpha.*...
        (n-j+1+alpha)).*f2(t(j-1),x(j-1),y(j-1),z(j-1)).*t(j-1).^(tau-1));
    
    z(n+1)=z(1)+(tau*(h^alpha)./gamma(alpha+2)).*sum(((n+1-j).^alpha.*...
        (n-j+2+alpha)-(n-j).^alpha.*(n-j+2+2*alpha)).*...
        f3(t(j),x(j),y(j),z(j)).*t(j).^(tau-1)-((n+1-j).^(alpha+1)-(n-j).^alpha.*...
        (n-j+1+alpha)).*f3(t(j-1),x(j-1),y(j-1),z(j-1)).*t(j-1).^(tau-1));

    t(n+1)=t(n)+h;

end


% Use following lines to simulate the phase portrait
curve=animatedline;
for k=1:length(t)
    addpoints(curve,y(k),z(k))
    drawnow
end