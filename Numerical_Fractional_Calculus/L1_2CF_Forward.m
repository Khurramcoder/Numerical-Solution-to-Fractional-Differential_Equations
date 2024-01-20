% Akman,T.,Yildiz,B.,& Baleanu,D.(2018).
% Numerical discretization of Caputo-Fabrizio Derivative.
% Computational & Applied Mathematics,37(3),3307-3333.


% Numerical L1-2 Method to approximate the Caputo Fabrizio-Forward Derivative

%Example f(x)=cos(4x) over [0,2].

clc;clear;close all;


%Inputs

h=0.1;
x(1)=0;
xlast=2;
x=x(1):h:xlast;
k=ceil((xlast-x(1))/h);


alpha =0.1; 

M=1;w=4;
lam=-alpha/(1-alpha);
f=@(x)cos(w*x);
% The Algorithm
j=2:k+1;


C=(1/(2*h))*(f((j-2)*h).*(1-2.*j))+(1/h).*(f((j-1)*h).*(2*j-2))+(1./(2*h)).*(f(j*h).*(3-2*j));
Ej=exp(lam*(k*h-(j-1)*h));
Ejm1=exp(lam*(k*h-(j-2)*h));
D=(1/h^2)*(f((j-2)*h)-2*f((j-1)*h)+f(j*h));
E=(j-1).*h.*Ej-((1-alpha)/alpha)*Ej-(j-2).*h.*Ejm1+((1-alpha)/alpha)*Ejm1;
L12CF_Forward=(1/alpha).*sum(C.*(Ej-Ejm1)+D.*E);


%Preparation for the Exact Solution

syms x

Exact=(M/(1-alpha))*int(diff(f,x)*exp(lam*(xlast-x)),x,0,2);

Exact=double(Exact);

Error = abs(Exact-L12CF_Forward);

disp('   Steps    Stepsize    Exact    Approximate    Error');


disp('-----------------------------------------------------');


Result = [k   h   Exact    L12CF_Forward      Error]