% Solis-Perez, J.E., Gomez-Aguilar, J.F., & Atangana, A. (2018). Novel
% numerical method for solving variable-order fractional differential
% equations with power, exponential and Mittag-Leffler laws. Chaos,
% Solitons & Fractal, 114, 175-185.

% This numerical scheme is for
% Atangana-Baleanu-Caputo(Mittag-Leffler)differential operator

clc; clear; close all;

%Inputs

h=0.01; t(1)=0; tfinal=200; t=t(1):h:tfinal; N=ceil((tfinal-t(1))/h);
% Initial Conditions
x(1)=2; y(1)=-1; z(1)=1;
% Constant fractional order
alpha=0.976;
%parameters of the model
a=1;b=0;c=1;

%Given below is the set of ODEs(Financial System given in the above paper
%by (3.1a-3.1c))
f1= @(t,x,y,z)z+(y-a).*x;
f2= @(t,x,y,z)1-b*y-x.^2;
f3= @(t,x,y,z)-x-c*z;

for n=1:N
    m=2:n;
    x(n+1)=x(1)+((gamma(alpha).*(1-alpha))/(gamma(alpha).*(1-alpha)+alpha)).*f1(t(n),x(n),y(n),z(n))+...
        ((h.^(alpha))/((alpha+1).*((1-alpha)).*(gamma(alpha))+alpha)).*...
        sum(((n+1-m).^alpha.*(n-m+2+alpha)-(n-m).^alpha.*(n-m+2+2.*alpha)).*f1(t(m),x(m),y(m),z(m))-...
        ((n+1-m).^(alpha+1)-(n-m).^alpha.*(n-m+1+alpha)).*f1(t(m-1),x(m-1),y(m-1),z(m-1)));
    
    y(n+1)=y(1)+((gamma(alpha).*(1-alpha))/(gamma(alpha).*(1-alpha)+alpha)).*f2(t(n),x(n),y(n),z(n))+...
        ((h.^(alpha))/((alpha+1).*((1-alpha)).*(gamma(alpha))+alpha)).*...
        sum(((n+1-m).^alpha.*(n-m+2+alpha)-(n-m).^alpha.*(n-m+2+2.*alpha)).*f2(t(m),x(m),y(m),z(m))-...
        ((n+1-m).^(alpha+1)-(n-m).^alpha.*(n-m+1+alpha)).*f2(t(m-1),x(m-1),y(m-1),z(m-1)));
    
    z(n+1)=z(1)+((gamma(alpha).*(1-alpha))/(gamma(alpha).*(1-alpha)+alpha)).*f3(t(n),x(n),y(n),z(n))+...
        ((h.^(alpha))/((alpha+1).*((1-alpha)).*(gamma(alpha))+alpha)).*...
        sum(((n+1-m).^alpha.*(n-m+2+alpha)-(n-m).^alpha.*(n-m+2+2.*alpha)).*f3(t(m),x(m),y(m),z(m))-...
        ((n+1-m).^(alpha+1)-(n-m).^alpha.*(n-m+1+alpha)).*f3(t(m-1),x(m-1),y(m-1),z(m-1)));
    
    t(n+1)=t(n)+h;
end
figure(1)
plot(x,y)
xlabel('x'),ylabel('y'),zlabel('z'),legend('Constant Fractional-order')


% Variable version of the Atangana-Baleanu-Caputo(Mittag-Leffler) Algorithm Starts


alpha = @(t) 0.97+0.03*cos(t/10);

for n=1:N
    m=2:n;
    x(n+1)=x(1)+((gamma(alpha(t(m))).*(1-alpha(t(m))))/(gamma(alpha(t(m))).*(1-alpha(t(m)))+alpha(t(m)))).*f1(t(n),x(n),y(n),z(n))+...
        ((h.^(alpha(t(m))))/((alpha(t(m))+1).*((1-alpha(t(m)))).*(gamma(alpha(t(m))))+alpha(t(m)))).*...
        sum(((n+1-m).^alpha(t(m)).*(n-m+2+alpha(t(m)))-(n-m).^alpha(t(m)).*(n-m+2+2.*alpha(t(m)))).*f1(t(m),x(m),y(m),z(m))-...
        ((n+1-m).^(alpha(t(m))+1)-(n-m).^alpha(t(m)).*(n-m+1+alpha(t(m)))).*f1(t(m-1),x(m-1),y(m-1),z(m-1)));
    
    y(n+1)=y(1)+((gamma(alpha(t(m))).*(1-alpha(t(m))))/(gamma(alpha(t(m))).*(1-alpha(t(m)))+alpha(t(m)))).*f2(t(n),x(n),y(n),z(n))+...
        ((h.^(alpha(t(m))))/((alpha(t(m))+1).*((1-alpha(t(m)))).*(gamma(alpha(t(m))))+alpha(t(m)))).*...
        sum(((n+1-m).^alpha(t(m)).*(n-m+2+alpha(t(m)))-(n-m).^alpha(t(m)).*(n-m+2+2.*alpha(t(m)))).*f2(t(m),x(m),y(m),z(m))-...
        ((n+1-m).^(alpha(t(m))+1)-(n-m).^alpha(t(m)).*(n-m+1+alpha(t(m)))).*f2(t(m-1),x(m-1),y(m-1),z(m-1)));
    
    z(n+1)=z(1)+((gamma(alpha(t(m))).*(1-alpha(t(m))))/(gamma(alpha(t(m))).*(1-alpha(t(m)))+alpha(t(m)))).*f3(t(n),x(n),y(n),z(n))+...
        ((h.^(alpha(t(m))))/((alpha(t(m))+1).*((1-alpha(t(m)))).*(gamma(alpha(t(m))))+alpha(t(m)))).*...
        sum(((n+1-m).^alpha(t(m)).*(n-m+2+alpha(t(m)))-(n-m).^alpha(t(m)).*(n-m+2+2.*alpha(t(m)))).*f3(t(m),x(m),y(m),z(m))-...
        ((n+1-m).^(alpha(t(m))+1)-(n-m).^alpha(t(m)).*(n-m+1+alpha(t(m)))).*f3(t(m-1),x(m-1),y(m-1),z(m-1)));
    t(n+1)=t(n)+h;
end
figure(2)
plot(x,y)
xlabel('x'),ylabel('y'),zlabel('z'),legend('Variable fractional-order')