%Jajarmi, A., & Baleanu, D.(2018). A new fractional analysis on the interaction of HIV
%with CD4+ T-cells, Chaos, Solitons & Fractals, 113,221-229.
%The numerical method coded in this m file is given in Eq(23)in the above
%research paper

clc;clear;close all;
%Inputs
h=0.01; t(1)=0; tfinal=200; t=t(1):h:tfinal; N=ceil((tfinal-t(1))/h);
%Initial Conditions
x(1)=2; y(1)=1; z(1)=1;
alpha=0.964;%fractional order
g=9.91;A=2.1;B=0.24;D=0;b=12.57;%Parameters
%The #D dynamical system
f1=@(t,x,y,z) g*(y-B*sin(pi*x/(2*A)+D));
f2=@(t,x,y,z)x-y+z;
f3=@(t,x,y,z)-b*y;
% Caputo-Fabrizio Euler Method starts
tic;
for n=1:N
j=1:n;
x(n+1)=x(1)+((1-alpha)*f1(t(n),x(n),y(n),z(n)))+(alpha.*h).*sum(f1(t(j),x(j),y(j),z(j)));
y(n+1)=y(1)+((1-alpha)*f2(t(n),x(n),y(n),z(n)))+(alpha.*h).*sum(f2(t(j),x(j),y(j),z(j)));
z(n+1)=z(1)+((1-alpha)*f3(t(n),x(n),y(n),z(n)))+(alpha.*h).*sum(f3(t(j),x(j),y(j),z(j)));
t(n+1)=t(n)+h;
end
toc;
figure(1)
plot(x,y),xlabel('x'),ylabel('y'),
legend('Phase Plane')
figure(2)
plot3(x,y,z),xlabel('x'),ylabel('y'),zlabel('z'),
legend('3D Surface')