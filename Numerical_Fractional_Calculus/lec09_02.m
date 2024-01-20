% Solis-Perez, J.E., Gomez-Aguilar, J.F., & Atangana, A. (2018). Novel
% numerical method for solving variable-order fractional differential
% equations with power, exponential and Mittag-Leffler laws. Chaos,
% Solitons & Fractal, 114, 175-185.


% This numerical scheme is for Caputo differential operator

clc; clear; close all;

%Inputs

h=0.1; t(1)=0; tfinal=300; t=t(1):h:tfinal; N=ceil((tfinal-t(1))/h);
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

% Constant version of the Caputo Algorithm Starts

tic;
for n=1:N
    j=2:n;
    x(n+1)=x(1)+((h.^alpha)./(alpha.*(1+alpha).*gamma(alpha))).*...
        sum(((n+1-j).^alpha.*(n-j+2+alpha)-(n-j).^alpha.*(n-j+2+2.*alpha)).*f1(t(j),x(j),y(j),z(j))-...
        ((n+1-j).^(alpha+1)-(n-j).^alpha.*(n-j+1+alpha)).*f1(t(j-1),x(j-1),y(j-1),z(j-1)));
    
    y(n+1)=y(1)+((h.^alpha)./(alpha.*(1+alpha).*gamma(alpha))).*...
        sum(((n+1-j).^alpha.*(n-j+2+alpha)-(n-j).^alpha.*(n-j+2+2.*alpha)).*f2(t(j),x(j),y(j),z(j))-...
        ((n+1-j).^(alpha+1)-(n-j).^alpha.*(n-j+1+alpha)).*f2(t(j-1),x(j-1),y(j-1),z(j-1)));
    
    z(n+1)=z(1)+((h.^alpha)./(alpha.*(1+alpha).*gamma(alpha))).*...
        sum(((n+1-j).^alpha.*(n-j+2+alpha)-(n-j).^alpha.*(n-j+2+2.*alpha)).*f3(t(j),x(j),y(j),z(j))-...
        ((n+1-j).^(alpha+1)-(n-j).^alpha.*(n-j+1+alpha)).*f3(t(j-1),x(j-1),y(j-1),z(j-1)));
    
    t(n+1) =t(n)+h;
end
toc;
figure(1)
plot(x,y)
xlabel('x'),ylabel('y'),zlabel('z'),legend('Constant Fractional-order')

alpha = @(t) 1./(1+exp(-t));


% Variable version of the Caputo Algorithm Starts

tic;
for n=1:N
    j=2:n;
    x(n+1)=x(1)+sum(((h.^(alpha(t(j))))./((alpha(t(j))).*(1+(alpha(t(j)))).*gamma((alpha(t(j)))))).*...
    ((n+1-j).^(alpha(t(j))).*(n-j+2+(alpha(t(j))))-(n-j).^(alpha(t(j))).*(n-j+2+2.*(alpha(t(j))))).*...
    f1(t(j),x(j),y(j),z(j))-((n+1-j).^((alpha(t(j)))+1)-(n-j).^(alpha(t(j))).*(n-j+1+(alpha(t(j))))).*...
    f1(t(j-1),x(j-1),y(j-1),z(j-1)).*((h.^(alpha(t(j))))./((alpha(t(j))).*(1+(alpha(t(j)))).*...
    gamma((alpha(t(j)))))));

    y(n+1)=y(1)+sum(((h.^(alpha(t(j))))./((alpha(t(j))).*(1+(alpha(t(j)))).*gamma((alpha(t(j)))))).*...
    ((n+1-j).^(alpha(t(j))).*(n-j+2+(alpha(t(j))))-(n-j).^(alpha(t(j))).*(n-j+2+2.*(alpha(t(j))))).*...
    f2(t(j),x(j),y(j),z(j))-((n+1-j).^((alpha(t(j)))+1)-(n-j).^(alpha(t(j))).*(n-j+1+(alpha(t(j))))).*...
    f2(t(j-1),x(j-1),y(j-1),z(j-1)).*((h.^(alpha(t(j))))./((alpha(t(j))).*(1+(alpha(t(j)))).*...
    gamma((alpha(t(j)))))));

    z(n+1)=z(1)+sum(((h.^(alpha(t(j))))./((alpha(t(j))).*(1+(alpha(t(j)))).*gamma((alpha(t(j)))))).*...
    ((n+1-j).^(alpha(t(j))).*(n-j+2+(alpha(t(j))))-(n-j).^(alpha(t(j))).*(n-j+2+2.*(alpha(t(j))))).*...
    f3(t(j),x(j),y(j),z(j))-((n+1-j).^((alpha(t(j)))+1)-(n-j).^(alpha(t(j))).*(n-j+1+(alpha(t(j))))).*...
    f3(t(j-1),x(j-1),y(j-1),z(j-1)).*((h.^(alpha(t(j))))./((alpha(t(j))).*(1+(alpha(t(j)))).*...
    gamma((alpha(t(j)))))));
t(n+1)=t(n)+h;

end
toc;
figure(2)
plot(x,y)
xlabel('x'),ylabel('y'),zlabel('z'),legend('Variable fractional-order')
