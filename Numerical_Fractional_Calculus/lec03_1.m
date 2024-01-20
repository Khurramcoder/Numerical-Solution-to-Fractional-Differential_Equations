% This numerical scheme is for Caputo Fractal-Fractional
% This scheme is proposed by Prof. Dr. Atangana.
% Code is for Model (3.4) in the research paper cited as:
% Atangana, A., & Qureshi, S. (2019). Modeling attractors of chaotic dynamical
%systems with fractal-fractional operators. Chaos, Solitons & Fractals, 123, 320-337.
clc; clear;close all;
h= 0.01; t(1)=0.1; tfinal=100; t-t(1):h:tfinal; N=ceil((tfinal-t(1))/h); % Inputs
x(1)=2; y(1)=1; z(1)=1; % Initial Conditions
alpha=0.97;% this is fractional order
tau=0.89; % this is fractal dimension
g = 8.91;A = 2.1;B = 0.24; D = 0; b=12.57;% Parameters in the Model(3.4)
% The Model (3.4) in the Research Paper
f1 = @(t,x,y,z) g*(y-B*sin(pi*x/(2*A)+D));
f2 = @(t,x,y,z) x-y+z;
f3= @(t,x,y,z) -b*y;
% Algorithm of the Caputo Fractal-Fractional starts
tic;
for n=1:N
j=2:n;
x(n+1)=x(1)+(tau*(h^alpha)./gamma(alpha+2)).*sum(((n+1-j).^alpha .*...
(n-j+2+alpha)-(n-j).^alpha.*(n-j+2+2*alpha)).*...
f1(t(j),x(j),y(j),z(j)).*t(j).^(tau-1)-((n+1-j).^(alpha+1)-(n-j).^alpha.*...
(n-j+1+alpha)).*f1(t(j-1),x(j-1),y(j-1),z(j-1)).*t(j-1).^(tau-1));
y(n+1)=y(1)+(tau*h^alpha./gamma(alpha+2)).*sum(((n+1-j).^alpha .*...
(n-j+2+alpha)-(n-j).^alpha .*(n-j+2+2*alpha)).*...
f2(t(j),x(j),y(j),z(j)).*t(j).^(tau-1)-((n+1-j).^(alpha+1)-(n-j).^alpha.*...
(n-j+1+alpha)).*f2(t(j-1),x(j-1),y(j-1),z(j-1)).*t(j-1).^(tau-1));
z(n+1)=z(1)+(tau*h^alpha./gamma(alpha+2)).*sum(((n+1-j).^alpha.*...
(n-j+2+alpha)-(n-j).^alpha .* (n-j+2+2*alpha)).*...
f3(t(j),x(j),y(j),z(j)).*t(j).^(tau-1)-((n+1-j).^(alpha+1)-(n-j).^alpha .*...
(n-j+1+alpha)).*f3(t(j-1),x(j-1),y(j-1),z(j-1)).*t(j-1).^(tau-1));
t(n+1)=t(n)+h;
end
toc;
plot(x,y),xlabel('x'),ylabel('y'), % Plot (a) of Figure 7 in the research paper
