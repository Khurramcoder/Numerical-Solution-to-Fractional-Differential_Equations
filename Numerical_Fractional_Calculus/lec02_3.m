% Numerical scheme proposed by Mekkaoui Toufik and Abdon Atangana
% This is the code for Example 1 from the following paper
% Toufik, M., & Atangana, A. (2017). New numerical approximation of fractional derivative with non-local and non-singular kernel:
%application to chaotic models. The European Physical Journal Plus, 132(10), 1-16.
clc;clear;close all;
format shorte
% Inputs
alpha=3/4;abc=1-alpha+alpha/gamma(alpha);
% Initial Conditions
t(1)=0; y(1)=0;h= 0.1; tfinal=1; t=t(1):h:tfinal;
% Number of Iterations
N=ceil((tfinal-t(1))/h);
% The given fractional-order ODE under the ABC operator
f = @(t,y) t.^2;
Exact = (1-alpha)/abc*t.^2+(2/(gamma(alpha)*abc*(alpha^2+3*alpha+2)))*t.^(alpha+2);
% ABC Algorithm starts
for n=1:N
k=2:n;
y(n+1)=y(1)+((1-alpha)/abc)*f(t(n),y(n))+(alpha/abc)*(h^alpha/gamma(alpha+2)).*...
sum(((n+1-k).^alpha.*(n-k+2+alpha)-(n-k).^alpha.*(n-k+2+2*alpha)).*f(t(k),y(k))-...
((n+1-k).^(alpha+1)-(n-k).^alpha.*(n-k+1+alpha)).*f(t(k-1),y(k-1)));
t(n+1)=t(n)+h;
end
Errors=abs(Exact-y);
Last_Error = Errors(end),
plot(t,Exact,'k-')
hold on
plot(t,y,'b-')