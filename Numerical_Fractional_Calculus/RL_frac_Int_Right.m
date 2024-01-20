% Diethelm,K.,& Karniadakis,G.E.(2019)
%Fundamental approaches for the numerical handling of fractional operators
%and time-fractional differential equations.
%Handbook of Fractional Calculus with Applications,3,1-22.

%Right Product Rectangle Rule to approximate Riemann Liouville Fractional
%Integral.

clc;clear;close all;
format compact
%Inputs

h=0.1;x(1)=0;xlast=1;x=x(1):h:xlast;

j=ceil((xlast-x(1))/h);

alpha=0.845;
f=@(x)sin(x);

k=2:j;
Weight=(h*alpha/gamma(alpha+1));
Right_RL_Int=Weight*sum(((j-k+1).^alpha-(j-k).^alpha).*f(x(k)));


%Exact Solution

syms q

Exact=eval(symsum((-1)^q/gamma(alpha+2*q+2),q,0,Inf));

Error=abs(Exact-Right_RL_Int);

disp('    Steps        Stepsize       Exact       Approximate    Error')
disp('-----------------------------------------------------------------')


Results=[   j       h          Exact Right_RL_Int              Error]