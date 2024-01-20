% Solis-Perez, J.E., Gomez-Aguilar, J.F., & Atangana, A. (2018). Novel
% numerical method for solving variable-order fractional differential
% equations with power, exponential and Mittag-Leffler laws. Chaos,
% Solitons & Fractal, 114, 175-185.



% This numerical scheme is for Caputo differential operator

clc; clear; close all;

%Inputs

h=0.01; t(1)=0; tfinal=1000; t=t(1):h:tfinal; N=ceil((tfinal-t(1))/h);
% Initial Conditions
x(1)=2; y(1)=-1; z(1)=1;
% Constant fractional order
alpha = 0.976;
% Parameters of the model
a=1; b=0.1; c=1;
% Given below is the set of ODEs(Financial System given in the above paper
% by (3.1a-3.1c)
f1= @(t,x,y,z)z+(y-a).*x;
f2= @(t,x,y,z)1-b*y-x.^2;
f3= @(t,x,y,z)-x-c*z;

x(2)=x(1)+((h^alpha)/(gamma(alpha+1))).*f1(t(1),x(1),y(1),z(1));
y(2)=y(1)+((h^alpha)/(gamma(alpha+1))).*f2(t(1),x(1),y(1),z(1));
z(2)=z(1)+((h^alpha)/(gamma(alpha+1))).*f3(t(1),x(1),y(1),z(1));

% Constant version of the Caputo-Fabrizio-Caputo Algorithm starts
for n=2:N
    x(n+1)=x(n)+(0.5*(2-alpha)*(1-alpha)+0.25*3*h*alpha*(2-alpha))*f1(t(n),x(n),y(n),z(n))-...
        (0.5*(2-alpha)*(1-alpha)+0.25*h*alpha*(2-alpha))*f1(t(n-1),x(n-1),y(n-1),z(n-1));
    
    y(n+1)=y(n)+(0.5*(2-alpha)*(1-alpha)+0.25*3*h*alpha*(2-alpha))*f2(t(n),x(n),y(n),z(n))-...
        (0.5*(2-alpha)*(1-alpha)+0.25*h*alpha*(2-alpha))*f2(t(n-1),x(n-1),y(n-1),z(n-1));
    
    z(n+1)=z(n)+(0.5*(2-alpha)*(1-alpha)+0.25*3*h*alpha*(2-alpha))*f3(t(n),x(n),y(n),z(n))-...
        (0.5*(2-alpha)*(1-alpha)+0.25*h*alpha*(2-alpha))*f3(t(n-1),x(n-1),y(n-1),z(n-1));


    t(n+1)=t(n)+h;
end

figure(1)
plot3(x,y,z)
xlabel('x'),ylabel('y'),zlabel('z'),legend('Constant Fractional-order')

alpha = @(t) 1./(1+exp(-t));

% Variable version of the Caputo-Fabrizio-Caputo Algorithm Starts
for n=2:N
    x(n+1)=x(n)+(0.5*(2-alpha(t(n)))*(1-alpha(t(n)))+0.25*3*h*alpha(t(n))*(2-alpha(t(n))))*...
        f1(t(n),x(n),y(n),z(n))-(0.5*(2-alpha(t(n)))*(1-alpha(t(n)))+0.25*h*alpha(t(n))*...
        (2-alpha(t(n))))*f1(t(n-1),x(n-1),y(n-1),z(n-1));
    
    y(n+1)=y(n)+(0.5*(2-alpha(t(n)))*(1-alpha(t(n)))+0.25*3*h*alpha(t(n))*(2-alpha(t(n))))*...
        f2(t(n),x(n),y(n),z(n))-(0.5*(2-alpha(t(n)))*(1-alpha(t(n)))+0.25*h*alpha(t(n))*...
        (2-alpha(t(n))))*f2(t(n-1),x(n-1),y(n-1),z(n-1));
    
    z(n+1)=z(n)+(0.5*(2-alpha(t(n)))*(1-alpha(t(n)))+0.25*3*h*alpha(t(n))*(2-alpha(t(n))))*...
        f3(t(n),x(n),y(n),z(n))-(0.5*(2-alpha(t(n)))*(1-alpha(t(n)))+0.25*h*alpha(t(n))*...
        (2-alpha(t(n))))*f3(t(n-1),x(n-1),y(n-1),z(n-1));
    
t(n+1)=t(n)+h;
end
figure(2)
plot3(x,y,z)
xlabel('x'),ylabel('y'),zlabel('z'),legend('Variable fractional-order')


