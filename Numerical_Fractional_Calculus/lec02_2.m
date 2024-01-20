% Numerical scheme proposed by Mekkaoui Toufik and Abdon Atangana
% This is the code for Example 3 from the following paper
% Toufik, M., & Atangana, A. (2017). New numerical approximation of fractional derivative with non-local and non-singular kernel:
%application to chaotic models. The European Physical Journal Plus, 132(10), 1-16.
clc;clear;close all;
% Inputs
h= 0.01; t(1)=0; x(1)=1; y(1)=1; z(1)=1;
a=0.998; abc=1-a+a/gamma(a); R=60;
tfinal=100; t=t(1):h:tfinal; N=ceil((tfinal-t(1))/h);
% Differential Equations
f1 = @(t,x,y,z) y-x;
f2 = @(t,x,y,z) -z.*tanh(x);
f3 = @(t,x,y,z) -R+x.*y+abs(y);
% ABC Algorithm starts
for n=1:N
k=2:n;
x(n+1)=x(1)+((1-a)/abc)*f1(t(n),x(n),y(n),z(n))+(a/abc)*(h^a/gamma(a+2)).*...
sum(((n+1-k).^a.*(n-k+2+a)-(n-k).^a .*(n-k+2+2*a)).*f1(t(k),x(k),y(k),z(k))-....
((n+1-k).^(a+1)-(n-k).^a.*(n-k+1+a)).*f1(t(k-1),x(k-1),y(k-1),z(k-1)));
y(n+1)=y(1)+((1-a)/abc)*f2(t(n),x(n),y(n),z(n))+(a/abc)*(h^a/gamma(a+2)).*...
sum(((n+1-k).^a.*(n-k+2+a)-(n-k).^a .*(n-k+2+2*a)).*f2(t(k),x(k),y(k),z(k))-...
((n+1-k).^(a+1)-(n-k).^a.*(n-k+1+a)).*f2(t(k-1),x(k-1),y(k-1),z(k-1)));
z(n+1)=z(1)+((1-a)/abc)*f3(t(n),x(n),y(n),z(n))+(a/abc)*(h^a/gamma(a+2)).*... 
sum(((n+1-k).^a.*(n-k+2+a)-(n-k).^a.*(n-k+2+2*a)).*f3(t(k),x(k),y(k),z(k))-....
((n+1-k).^(a+1)-(n-k).^a .* (n-k+1+a)).*f3(t(k-1),x(k-1),y(k-1),z(k-1)));
t(n+1)=t(n)+h;
end
% Phase portrait
plot(x,z)
xlabel('x')
ylabel('z')
