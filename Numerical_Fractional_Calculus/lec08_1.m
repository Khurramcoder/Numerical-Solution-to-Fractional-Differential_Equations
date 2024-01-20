%Numerical Caputo Fractional Integral
%The Caputo Fractional Derivative Eq(3.4)
%The Research Paper is as follows:
%Odibat, Z.(2006). Approximations of fractional integrals and Caputo fractional derivatives.
%Application Mathematics and Computations, 178(2), 527-533.

%Example 2 f(x)=sin(x) over [0,1].


clc; clear; close all;

%Inputs

h=0.1; x(1)=0; xlast=1; x=x(1):h:xlast; k=ceil((xlast-x(1))/h);
alpha=0.5;
% The Algorithm
tic;
j=1:k-1;
Caputo_D = ((h^(1-alpha))/gamma(3-alpha)).*(((k-1)^(2-alpha)-(k+alpha-2)*k^(1-alpha))*...
    cos(x(1))+cos(xlast)+sum(((k-j+1).^(2-alpha)-2*(k-j).^(2-alpha)+(k-j-1).^...
    (2-alpha)).*cos(x(j+1))));
toc;

%Exact Solution
syms i
actual=eval(symsum((-1)^i/gamma(2*i+2-0.5),i,0,Inf));
% Finding Absolute Error
Error =abs(actual-Caputo_D);
disp('   Steps          Stepsize   Approximate   Error')
disp('-------------------------------------------------')
Results=[k h Caputo_D Error]




