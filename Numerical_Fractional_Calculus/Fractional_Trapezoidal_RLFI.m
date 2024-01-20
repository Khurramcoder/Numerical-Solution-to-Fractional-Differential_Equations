%Dr. Sania Qureshi


% Fractional Trapezoidal Rule to approximate Riemann Liouville Fractional
% Integral

clc;clear;close all;
format compact
%Inputs

h=0.1;x(1)=0;xlast=1;x=x(1):h:xlast;

n=ceil((xlast-x(1))/h);

alpha=0.95;

f=@(x)sin(x);


% The Algorithm

k=2:n-1;

First=(n-1)^(alpha+1)-(n-1-alpha)*(n^alpha);

Middle=(n-k+1).^(alpha+1)+(n-1-k).^(alpha+1)-2.*(n-k).^(alpha+1);

Coeff=h^alpha/gamma(alpha+2);

Trap_RL_Int=Coeff*(First*f(x(1))+sum((Middle.*f(x(k))))+f(xlast));


% Exact Solution

syms q


Exact=eval(symsum((-1)^q/gamma(alpha+2*q+2),q,0,Inf));

Error=abs(Exact-Trap_RL_Int);

disp('    Steps        Stepsize       Exact       Approximate    Error')
disp('-----------------------------------------------------------------')


Results=[   n       h          Exact   Trap_RL_Int              Error]