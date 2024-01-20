% The Fractional Forward Euler's Method (22.05.2022)
clc; clear; close all;
% The inputs,
h=1/10;t(1)=0; y(1)=0; alpha=0.1; tfinal=1;t=t(1):h:tfinal;N=ceil((tfinal-t(1))/h);
% Exact Solution
Exact=t.^4.*mlf(alpha, 5,-(t.^alpha));
%Fractional-Order ODE
f=@(t,y) -y+(1/gamma(5-alpha)).*t.^(4-alpha);
% The Fractional Forward Euler Method
for n=1:N
j=1:n;
t(n+1)=t(n)+h;
y(n+1)=y(1)+((h^alpha)/(gamma(alpha+1))).*sum(((n-j+1).^(alpha)-(n-j).^(alpha)).*f(t(j),y(j)));
end
%Absolute Errors
Errors=abs(Exact-y);
Last_Error=Errors(end),