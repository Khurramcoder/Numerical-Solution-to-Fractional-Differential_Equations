% Diethelm,K.,& Karniadakis,G.E.(2019)
%Fundamental approaches for the numerical handling of fractional operators
%and time-fractional differential equations.
%Handbook of Fractional Calculus with Applications,3,1-22.

%Left Product Rectangle Rule to approximate Riemann Liouville Fractional
%Integral.

clc;clear;close all;

%Inputs

h=0.1;x(1)=0;xlast=1;x=x(1):h:xlast;

j=ceil((xlast-x(1))/h);

alpha=0.795;
f=@(x)sin(x);

k=1:j-1;
Weight=(h*alpha/gamma(alpha+1));
Left_RL_Int=Weight*sum(((j-k).*alpha-(j-k-1).^alpha).*f(x(k)));


%Exact Solution

syms q

Exact=eval(symsum((-1)^q/gamma(alpha+2*q+2),q,0,Inf));

Error=abs(Exact-Left_RL_Int);

disp('    Steps        Stepsize       Exact       Approximate    Error')
disp('-----------------------------------------------------------------')


Results=[j h Exact Left_RL_Int Error]

