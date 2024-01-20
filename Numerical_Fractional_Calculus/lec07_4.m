%Numerical Riemann Liouville Fractional Ibtegral
%Modified Trapezoidal Rule Eq(2.1)
%The Research Paper is as follows:
%Odibat, Z.(2006). Approximations of fractional integrals and Caputo fractional derivatives.
%Application Mathematics and Computations, 178(2), 527-533.


%Example 1 f(x)=sin(x) over [0,1].


clc; clear; close all;

%Inputs

h=0.1; x(1)=0; xlast=1; x=x(1):h:xlast; k=ceil((xlast-x(1))/h);
alpha=0.5;
f0=sin(x(1));fa=sin(xlast);

%The Algorithm

tic;

j=1:k-1;

RL_Int= (h^alpha/gamma(alpha+2))*(((k-1)^(alpha+1)-(k-alpha-1)*k^alpha)*f0+...
    fa+sum(((k-j+1).^(alpha+1)-2*(k-j).^(alpha+1)+(k-j-1).^...
    (alpha+1)).*sin(x(j+1))));
toc;

%Exact Solution
syms i
ExactI=eval(symsum((-1)^i/gamma(alpha+2*i+2),i,0,Inf));
Error =abs(ExactI-RL_Int);
disp('   Steps          Stepsize   Approximate   Error')
disp('-------------------------------------------------')
Results=[k h RL_Int Error]

